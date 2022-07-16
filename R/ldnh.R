## lndh: an R package for processing square-wave voltammetry data
## Created by: Stephen Ramsey, Oregon State University
## License and copyright:  see file `LICENSE` and see
## documentation in github:ramseylab/ldnh-private

## ---------- GLOBAL DEFAULT PARAMETER VALUES ---------
detilt_start_v_default <- 1.00
detilt_end_v_default <- 1.15

peakfind_start_v_default <- 1.00
peakfind_end_v_default <- 1.10
## ----------------------------------------------------

#'
#' @import xlsx remotes Rdistance magrittr ggplot2 SplinesUtils  
#'
function() {}  ## this is required in order to force roxygen2::roxygenise to
               ## generate "import" directives in the NAMESPACE file


plot_conc_faceted_voltammograms <- function(df, cur_var, ylab, file_name) {
    min_voltage <- min(df$potential)
    max_voltage <- max(df$potential)
    p <- ggplot2::ggplot(data=df,
                         ggplot2::aes_string(x="potential", y=cur_var, colour="device")) +
        ggplot2::facet_grid(ggplot2::vars(conc_factor)) +
        ggplot2::ylab(ylab) +
        ggplot2::xlim(min_voltage, max_voltage) +
        ggplot2::theme_classic() +
        ggplot2::geom_line(size=0.5)
    ggplot2::ggsave(p, file=file_name, height=6, width=3)
}

make_residual_calculator_passthrough <- function(current_var) {
    function(df) {
        df$resid <- df[[current_var]]
        df
    }
}

## convert each element of a vector of "signal" values (each value corresponding
## to a single voltammogram) to an estimated carbamazepine concentration, using
## a linear calibration curve fit (carbamazepine concentration 
transform_signals_to_conc <- function(signal_values, model_coefficients) {
    fit_slope <- model_coefficients["signal"]
    fit_intercept <- model_coefficients["(Intercept)"]
    conc_values <- fit_intercept + signal_values*fit_slope
}

detilt_linear_model <- function(df, input_y_variable_name, output_y_variable_name, potential_range) {
    start_potential <- potential_range[1]
    end_potential <- potential_range[2]
    stopifnot(end_potential > start_potential)

    lapply(levels(df$conc_factor), function(level) {
        dat_sub <- subset(df, conc_factor==level)
        lapply(levels(dat_sub$device), function(dev) {
            dat_sub_sub <- subset(dat_sub, device==dev) 
            dat_sub_sub_censored <- subset(dat_sub_sub, potential <= start_potential |
                                                        potential >= end_potential)
            spline_model <- SplinesUtils::SmoothSplineAsPiecePoly(smooth.spline(dat_sub_sub_censored$potential,
                                                                                dat_sub_sub_censored[[input_y_variable_name]]))
            dat_sub_sub[[output_y_variable_name]] <- dat_sub_sub[[input_y_variable_name]] - predict(spline_model, dat_sub_sub$potential)
            dat_sub_sub
        }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .) -> df

    df
}

is_unix <- function() {
    .Platform$OS.type == "unix"
}

find_stationary_points <- function(x, y) {
    ## from SO:48720040
    sm <- smooth.spline(x, y)
    oo <- SplinesUtils::SmoothSplineAsPiecePoly(sm)
    {
        if (is_unix()) {
            sink("/dev/null")
        } else {
            ## we must be on windows
            sink("nul")
        }
        xs <- solve(oo, deriv=1)
        sink()
    }
    ys <- predict(oo, xs)
    list(x=xs, y=ys, oo=oo, sm=sm)
}


find_and_mark_peaks <- function(df, input_y_var_name, peak_x_var_name, window_region) {
    window_left <- window_region[1]
    window_right <- window_region[2]
    window_center <- 0.5*(window_left + window_right)
       
    ## adjust the voltammograms so that they have zero signal value at the "window center"
    lapply(levels(df$conc_factor), function(level) {
        dat_sub <- subset(df, conc_factor==level)
        lapply(levels(dat_sub$device), function(dev) {
            dat_sub_sub <- subset(dat_sub, device==dev)
            sp_list <- find_stationary_points(dat_sub_sub$potential,
                                              dat_sub_sub[[input_y_var_name]])

            x_use <- NA
            
            if (length(sp_list$x) == 0) {
                cat(sprintf("no stationary points found in the whole potential range; conc=%0.1f; device=%d\n",
                              dat_sub_sub$conc[1],
                            dat_sub_sub$device[1]))
            } else {
                sp_list_subset <- subset(data.frame(sp_list[c("x", "y")]),
                                         x >= window_left &
                                         x <= window_right)
                if (nrow(sp_list_subset) == 1) {
                    x_use <- sp_list_subset$x
                    cat(sprintf("one stationary point within the default window; conc=%0.1f; device=%d; x=%0.3f\n",
                                dat_sub_sub$conc[1],
                                dat_sub_sub$device[1],
                                x_use))
                } else {
                    if (nrow(sp_list_subset) == 0) {
                        cat(sprintf("no stationary points found in the default window; conc=%0.1f; device=%d\n",
                                    dat_sub_sub$conc[1],
                                    dat_sub_sub$device[1]))
                    } else {
                        sec_deriv_func <- function(x) {
                            predict(sp_list$sm, x)$y
                        }
                        hvals <- sapply(sp_list_subset$x, function(x) {
                            Rdistance::secondDeriv(x, FUN=sec_deriv_func)})
                        min_ind <- which.min(hvals)
                        y_use <- sp_list_subset$y[min_ind]
                        x_use <- sp_list_subset$x[min_ind]
                        cat(sprintf("out of %d stationary points, picking the one with most negative 2nd deriv; conc=%0.1f; device=%d; at x=%0.3f\n",
                                    length(hvals),
                                    dat_sub_sub$conc[1],
                                    dat_sub_sub$device[1],
                                    x_use))
                    }
                }
            }

            dat_sub_sub$peak_potential <- x_use
            dat_sub_sub
        }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .) -> df

    avg_peak_center <- mean(df$peak_potential, na.rm=TRUE)
    if (is.nan(avg_peak_center)) {
        cat(sprintf("No samples had a stationary point within the window; using the window center %0.3f\n", window_center))
        avg_peak_center <- window_center
    }

    lapply(levels(df$conc_factor), function(level) {
        dat_sub <- subset(df, conc_factor==level)
        lapply(levels(dat_sub$device), function(dev) {
            dat_sub_sub <- subset(dat_sub, device==dev)
            dat_sub_sub[[peak_x_var_name]] <- if (! is.na(dat_sub_sub$peak_potential[1])) {
                                                    dat_sub_sub$peak_potential[1]
                                                } else {
                                                    avg_peak_center
                                                }
            dat_sub_sub
        }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
}

make_signal_extractor_hessian <- function(eval_hessian_at_potential_var_name,
                                          input_y_var_name,
                                          output_y_var_name) {
    function(df) {
        sub_df <- subset(df, in_window)
        sm <- smooth.spline(sub_df$potential, sub_df[[input_y_var_name]])
        window_center <- sub_df[[eval_hessian_at_potential_var_name]][1]
        res_df <- df[1, c("conc_factor", "device", "conc"), drop=FALSE]
        res_df[["peak_potential"]] <- sub_df[["peak_potential"]][1]
        res_df[[output_y_var_name]] <- as.vector(-Rdistance::secondDeriv(window_center, FUN=function(x) {predict(sm, x)$y}))
        res_df
    }
}

window_data <- function(df, peak_x_variable_name, new_variable_name, window_width) {
    lapply(levels(df$conc_factor), function(level) {
        dat_sub <- subset(df, conc_factor==level)
        lapply(levels(dat_sub$device), function(dev) {
            dat_sub_sub <- subset(dat_sub, device==dev)
            window_center <- dat_sub_sub[1, peak_x_variable_name]
            window_left <- window_center - 0.5*window_width
            window_right <- window_center + 0.5*window_width
            dat_sub_sub[[new_variable_name]] <- FALSE
            within(dat_sub_sub, {
                temp_var <- get(new_variable_name)
                temp_var[potential >= window_left & potential <= window_right] <- TRUE
                assign(new_variable_name, temp_var)
                temp_var <- NULL})
        }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
}


fit_and_evaluate_calibration_curve <- function(results_list,
                                               plot_file_prefix=NULL) {

    res_with_signals <- results_list$all_results
    
    ## make a calibration line from all the data
    model_fit <- lm(conc~signal, data=res_with_signals)

    ## extract the fit coefficients from the calibration line
    model_coefficients <- coefficients(model_fit)

    ## estimate conc levels from the signals, using the calibration line
    conc_values <- transform_signals_to_conc(res_with_signals$signal, model_coefficients)

    ## compute prediction errors
    conc_errors <- conc_values - res_with_signals$conc

    ## need to handle the zero-analyte samples specially (use intercept value for those)
    nz <- which(res_with_signals$conc_factor != "0")

    ## compute the relative prediction errors
    conc_errors_rel <- c(conc_errors[! nz] / abs(model_coefficients["(Intercept)"]),
                        conc_errors[nz] / res_with_signals$conc[nz])

    ## compute the average relative prediction errors, using L1 and L2 norms
    pct_error_avg_l1 <- sqrt(mean(abs(conc_errors_rel)))
    pct_error_avg_l2 <- sqrt(mean(conc_errors_rel * conc_errors_rel))

    ## compute the correlation coefficient
    corrcoef <- cor(res_with_signals$conc, res_with_signals$signal)
    
    ## make plots if a plotting filename prefix is provided
    if (! is.null(plot_file_prefix)) {
        p <- ggplot2::ggplot(data=res_with_signals,
                             ggplot2::aes(conc, signal, colour=device)) +
            ggplot2::theme_classic() +
            ggplot2::geom_point()
        ggplot2::ggsave(p, file=paste(plot_file_prefix, "-dot-plot.pdf", sep=""), width=4, height=3)        
    }

    c(results_list,
        list(avg_rel_err_l1=pct_error_avg_l1,
             avg_rel_err_l2=pct_error_avg_l2,
             r=corrcoef,
             r2=corrcoef*corrcoef,
             slope=model_coefficients["signal"],
             intercept=model_coefficients["(Intercept)"]))
}


## For a melted data frame "raw_df_fac" of voltammograms for various analyte
## concentrations, and given various hyperparameter values, compute a "signal"
## value from the residuals produced by the `residuals calculator` function
## passed as an arg; find the best-fit linear model relating the analyte
## concentration to signal values; compute the Pearson product-moment
## coefficient between the signal values and the annotated carbamzepine levels,
## and also back-convert each signal value to an estimated analyte level and
## compute the average normalized error between the estimated and actual analyte
## concentration, across all analyte concentrations and replicates; return the
## results as a list. 
analyze_fitter_and_extractor <- function(raw_df_fac,
                                         make_signal_extractor,
                                         signal_extractor_params,
                                         curr_var,
                                         plot_file_prefix=NULL) {

    ## make the function that extracts the signal
    signal_extractor <- do.call(make_signal_extractor,
                                signal_extractor_params)

    raw_fac_with_resid <- raw_df_fac %>%
        split(list(conc_factor=raw_df_fac$conc_factor, device=raw_df_fac$device)) 

    ## stack the data frames (one for each voltammogram) into a single tall data frame
    raw_fac_with_resid_df <- do.call(rbind, raw_fac_with_resid)

    ## make a data frame with the conc level, the device number, and the final signal value
    res_with_signals <- raw_fac_with_resid %>%
        lapply(signal_extractor) %>%
        do.call(rbind, .)

    if (! is.null(plot_file_prefix)) {
       plot_conc_faceted_voltammograms(raw_fac_with_resid_df, "neg_current", "signal",
                                       paste(plot_file_prefix,
                                             "-voltammograms.pdf",
                                             sep=""))

       plot_conc_faceted_voltammograms(raw_fac_with_resid_df, "log_neg_current", "signal",
                                       paste(plot_file_prefix,
                                             "-voltammograms-log.pdf",
                                             sep=""))
         
       plot_conc_faceted_voltammograms(raw_fac_with_resid_df, "rel_log_neg_current", "signal",
                                       paste(plot_file_prefix,
                                             "-voltammograms-log-detilted.pdf",
                                             sep=""))
    }
    
    list(all_results=res_with_signals)
}

#' Run the script `process-voltammograms.R` that is installed in the `ldnh` package.
#'
#' @param fit_calibration A logical indicating whether or not a calibration
#'     curve should be fit to the signal-vs-concentration data; if TRUE, then
#'     the `-processed.xlsx` output spreadsheet is generated; if FALSE, then the
#'     `-processed.xlsx` output spreadsheet is not generated. Default: TRUE
#' @param file_name An optional string specifying the location of the
#'     spreadsheet (xlsx format) of melted data to be processed.
#' @export
#' @examples
#' process_voltammograms()
#' process_voltammograms(fit_calibration=FALSE)
#' process_voltammograms(fit_calibration=TRUE, file_name="melted-data-20220429.xlsx")

process_voltammograms <- function(fit_calibration=TRUE, file_name=NULL) {

    args <- commandArgs(trailingOnly=TRUE)
    if (length(args) > 0 || ! is.null(file_name)) {
        file_spec <- if (length(args) > 0) {
                         args[1]
                     } else {
                         file_name
                     }
        cat(sprintf("Loading input from file: %s\n", file_spec))
    } else {

        cat("Specify the file (xlsx format) of melted data to load: ")
        file_spec <- readLines(n=1)

        if (nchar(file_spec)==0) {
            stop("response to prompt was empty; expected a nonzero-length response")
        }
    }

    if (! file.exists(file_spec)) {
        stop(sprintf("File does not exist: %s", file_spec))
    }

    cat("Load worksheet number: [1] ")

    sheet_index <- readLines(n=1)

    if (nchar(sheet_index)==0) {
        sheet_index <- 1
    } else {
        sheet_index <- as.integer(sheet_index)
    }

    cat(sprintf("Specify the starting voltage for the detilting procedure, as a floating-point number [default: %f] ", detilt_start_v_default))
    detilt_start_v <- readLines(n=1)
    if (nchar(detilt_start_v) == 0) {
        detilt_start_v <- detilt_start_v_default
    } else {
        detilt_start_v <- as.numeric(detilt_start_v)
    }

    cat(sprintf("Specify the ending voltage for the detilting procedure, as a floating-point number [default: %f] ", detilt_end_v_default))
    detilt_end_v <- readLines(n=1)
    if (nchar(detilt_end_v) == 0) {
        detilt_end_v <- detilt_end_v_default
    } else {
        detilt_end_v <- as.numeric(detilt_end_v)
    }


    cat(sprintf("Specify the starting voltage for the peakfinding procedure, as a floating-point number [default: %f] ", peakfind_start_v_default))
    peakfind_start_v <- readLines(n=1)
    if (nchar(peakfind_start_v) == 0) {
        peakfind_start_v <- peakfind_start_v_default
    } else {
        peakfind_start_v <- as.numeric(peakfind_start_v)
    }

    cat(sprintf("Specify the ending voltage for the peakfinding procedure, as a floating-point number [default: %f] ", peakfind_end_v_default))
    peakfind_end_v <- readLines(n=1)
    if (nchar(peakfind_end_v) == 0) {
        peakfind_end_v <- peakfind_end_v_default
    } else {
        peakfind_end_v <- as.numeric(peakfind_end_v)
    }

    allowed_column_names <- c("potential",
                              "device",
                              "conc",
                              "current")

    {   ## just red the first line of the file, and check that the column names are correct
        df <- xlsx::read.xlsx(file_spec, sheetIndex=sheet_index, header=TRUE, endRow=2)

        if (! all(names(df) %in% allowed_column_names)) {
            stop(sprintf("Column name(s) in the file that are not recognized: %s",
                         paste(names(df)[which(! (names(df) %in% allowed_column_names))])))
        }
    }

    output_file_spec_prefix_default <- strsplit(file_spec, "\\.xlsx")[[1]]
    cat(sprintf("Specify the filespec prefix to be used for output files: [%s] ", output_file_spec_prefix_default))
    output_file_spec_prefix <- readLines(n=1)
    if (nchar(output_file_spec_prefix) == 0) {
        output_file_spec_prefix <- output_file_spec_prefix_default
    }

    log_file_name <- paste(output_file_spec_prefix, "-run.log", sep="")
    sink(log_file_name)
    on.exit(sink(NULL))

    writeLines(c("ldnh package: running process_voltammograms",
                 sprintf("Run started at: %s", as.character(Sys.time())),
                 as.character(sessionInfo()),
                 commandArgs(),
                 sprintf("Potential range to be used for detilting: %f - %f V", detilt_start_v, detilt_end_v),
                 sprintf("Potential range to be used for peakfinding: %f - %f V", peakfind_start_v, peakfind_end_v),
                 sprintf("User requested that a calibration curve be fit: %s", as.character(fit_calibration)),
                 sprintf("Input file: %s", file_spec)), sep="\n")

    xlsx::read.xlsx(file_spec, sheetIndex=sheet_index, header=TRUE) %>%  ## need magrittr pipe here
        within({
            ## turn the carbamazapine levels into factors rather than numeric types
            conc_factor <- factor(conc, levels=sort(unique(.$conc)))
            ## turn the replicate number values into factors instead of integer types
            device <- factor(device, levels=sort(unique(.$device)))
            ## create a column that is the negated current
            neg_current <- -current
            ## log2-transform the voltammograms
            log_neg_current <- log2(neg_current)
        }) |>
        
        ## de-tilt the voltammograms by subtracting the spline-fitted baseline
        detilt_linear_model("log_neg_current",
                            "rel_log_neg_current",
                            c(detilt_start_v, detilt_end_v)) |>

    ## adjust the voltammograms so that they have zero signal value at the "window center";
    ## we are going to search for a local maximum in the range [1.0, 1.1] volts
    find_and_mark_peaks("rel_log_neg_current",
                        "eval_hessian_at_potential",
                        c(peakfind_start_v, peakfind_end_v)) |>
    ## subset the data to a window of 0.1 volts width around the center that we found
    window_data("eval_hessian_at_potential", "in_window", 0.5*(peakfind_end_v - peakfind_start_v)) |>
    
    analyze_fitter_and_extractor(
        make_signal_extractor=make_signal_extractor_hessian,
        signal_extractor_params=list(eval_hessian_at_potential_var_name="eval_hessian_at_potential",
                                     input_y_var_name="rel_log_neg_current",
                                     output_y_var_name="signal"),
        curr_var="rel_log_neg_current_sub",
        plot_file_prefix=output_file_spec_prefix) -> results

    output_file_processed <- paste(output_file_spec_prefix, "-processed.xlsx", sep="")
    writeLines(c(sprintf("Saving processed data to file: %s", output_file_processed)), sep="\n")

    xlsx::write.xlsx(within(results$all_results, {conc_factor <- NULL}),
                     file=output_file_processed,
                     col.names=TRUE,
                     row.names=FALSE)

    if (fit_calibration) {

        results <- fit_and_evaluate_calibration_curve(results, 
                                                      plot_file_prefix=output_file_spec_prefix)
 
        output_file_summary <- paste(output_file_spec_prefix, "-summary.xlsx", sep="")
        writeLines(c(sprintf("Saving summary data to file: %s", output_file_summary)), sep="\n")
    
        xlsx::write.xlsx(data.frame(results[names(results)[which(names(results) != "all_results")]]),
                         file=output_file_summary,
                         col.names=TRUE,
                         row.names=FALSE)
    }

    writeLines(c(sprintf("Run completed at: %s", as.character(Sys.time()))), sep="\n")
}


get_voltammograms <- function(file_name) {
    read.table(file_name,
               sep=",",
               header=TRUE,
               skip=22) %>%
        setNames(c("potential", "current", "forward", "reverse")) %>%
        `[`(c("potential", "current"))
}

get_data <- function(file_name) {
    match_res <- stringr::str_match(file_name, "^(\\d+)_(\\d+)_(\\d+)_cbz([\\dp]+)_(\\d+).txt$")
    if (is.na(match_res[1, 1])) {
        stop(sprintf("Unable to process metadata from filename: %s", file_name))
    }
    year_str <- match_res[1, 2]
    if (nchar(year_str) == 2) {
        year_str <- paste("20", year_str, sep="")
    }
    year <- as.integer(year_str)
    month_str <- match_res[1, 3]
    month <- as.integer(month_str)
    day_str <- match_res[1, 4]
    day <- as.integer(day_str)
    melted_file_name <- sprintf("melted-data-%04d%02d%02d.xlsx", year, month, day)
    conc <- as.numeric(gsub("p", ".", match_res[, 5]))
    device_str <- match_res[, 6]
    device <- as.integer(device_str)
    metadata_df <- data.frame(melted_file_name=melted_file_name, device=device, conc=conc, row.names=file_name)
    v_df <- get_voltammograms(file_name)
    metadata_df %>% 
        replicate(n=nrow(v_df), ., simplify=FALSE) %>%
        do.call(rbind, .) %>%
        cbind(v_df) %>%
        within({current <- 1e+6 * current}) -> final_df

    final_df
}


#' Convert raw text data file from the potentiostat (one file per voltammogram)
#' to a single spreadsheet containing all of the voltammograms in a melted
#' form, with four columns: potential, device, conc, current. Current is in
#' microamperes.
#' @export
#' @examples
#' convert_raw_to_melted_xlsx()

convert_raw_to_melted_xlsx <- function() {
    all_files <- list.files(pattern=".txt$")
    writeLines(c("Files found: ", all_files),
               sep="\n")
    cat("Proceed? (answer y or n) [y]: ")
    response <- readLines(n=1)
    if (nchar(response) != 0 &&
        tolower(response) != "y") {
        cat("Exiting\n")
        return
    }

    all_files %>%
        lapply(get_data) %>%
        do.call(rbind, .) ->
        final_df

    stopifnot(! is.na(final_df$conc))
    stopifnot(! is.na(final_df$device))

    uniq_melted_file_name <- unique(final_df$melted_file_name)
    stopifnot(length(uniq_melted_file_name)==1)

    melted_file_name <- uniq_melted_file_name[1]

    final_df$melted_file_name <- NULL
    
    xlsx::write.xlsx(final_df[, c("potential", "device", "conc", "current")],
                     file=melted_file_name,
                     col.names=TRUE,
                     row.names=FALSE)

    cat(sprintf("Saved output as melted data spreadsheet: %s\n", melted_file_name))
}
