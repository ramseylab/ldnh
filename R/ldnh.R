#'
#' @import xlsx remotes Rdistance magrittr ggplot2 fork SplinesUtils
#'
function() {}  ## this is required in order to force roxygen2::roxygenise to
               ## generate "import" directives in the NAMESPACE file


detilt_start_v_default <- 0.5
detilt_end_v_default <- 0.9

peakfind_start_v_default <- 1.0
peakfind_end_v_default <- 1.1



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


detilt_linear_model <- function(df, new_variable_name, potential_range) {
    start_potential <- potential_range[1]
    end_potential <- potential_range[2]
    stopifnot(end_potential > start_potential)
    lmodel <- lm(log_neg_current ~ potential,
                data=subset(df, potential >= start_potential & potential <= end_potential))
    
    ## obtain the voltammograms of the zero-concentration (i.e., "blank") samples;
    avg_zero_conc_voltammogram <- subset(df, conc_factor=="0") %>%
    ## average the zero-conc voltammograms into a single representative voltammogram
        aggregate(cbind(current, neg_current, log_neg_current) ~ potential,
                  data=.,
                  FUN=mean)

    ## normalize log2 voltammograms to the average log2 voltammogram of the blank
    ## (i.e., zero conc) voltammogram
    df <- lapply(levels(df$conc_factor), function(level) {
        dat_sub <- subset(df, conc_factor==level)
        lapply(levels(dat_sub$device), function(dev) {
            dat_sub_sub <- subset(dat_sub, device==dev)
            dat_sub_sub[[new_variable_name]] <- dat_sub_sub$log_neg_current -
                avg_zero_conc_voltammogram$log_neg_current
            dat_sub_sub
        }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)

    df
}

find_stationary_points <- function(x, y) {
    ## from SO:48720040
    sm <- smooth.spline(x, y)
    oo <- SplinesUtils::SmoothSplineAsPiecePoly(sm)
    {
        sink("/dev/null")
        xs <- solve(oo, deriv=1)
        sink()
    }
    ys <- predict(oo, xs)
    list(x=xs, y=ys, oo=oo, sm=sm)
}


make_curves_kiss_avg <- function(df, new_variable_name, window_region) {
    window_left <- window_region[1]
    window_right <- window_region[2]
       
    ## adjust the voltammograms so that they have zero signal value at the "window center"
    lapply(levels(df$conc_factor), function(level) {
        dat_sub <- subset(df, conc_factor==level)
        lapply(levels(dat_sub$device), function(dev) {
            dat_sub_sub <- subset(dat_sub, device==dev)
            sp_list <- find_stationary_points(dat_sub_sub$potential,
                                              dat_sub_sub$rel_log_neg_current)
            use_default_window_center <- FALSE
            if (length(sp_list$x) == 0) {
                use_default_window_center <- TRUE
                cat(sprintf("no stationary points found in the whole potential range; conc=%0.1f; device=%d\n",
                              dat_sub_sub$conc[1],
                              dat_sub_sub$device[1]))
            } else {
                sp_list_subset <- subset(data.frame(sp_list[c("x", "y")]),
                                         x >= window_left &
                                         x <= window_right)
                if (nrow(sp_list_subset) == 1) {
                    y_use <- sp_list_subset$y
                    x_use <- sp_list_subset$x
                    cat(sprintf("one stationary point within the default window; conc=%0.1f; device=%d; y=%0.3f\n",
                                  dat_sub_sub$conc[1],
                                  dat_sub_sub$device[1],
                                  y_use))
                } else {
                    if (nrow(sp_list_subset) == 0) {
                        use_default_window_center <- TRUE
                        cat(sprintf("no stationary points found in the default window; conc=%0.1f; device=%d\n",
                                    dat_sub_sub$conc[1],
                                    dat_sub_sub$device[1]))
                    } else {
                        sec_deriv_func <- function(x) {
                            predict(sp_list$sm, x)$y
                        }
                        hvals <- sapply(sp_list_subset$x, function(x) {
                            Rdistance::secondDeriv(x, FUN=sec_deriv_func)})
                        min_ind <- which.max(abs(hvals))
                        y_use <- sp_list_subset$y[min_ind]
                        x_use <- sp_list_subset$x[min_ind]
                    }
                }
            }

            if (use_default_window_center) {
                x_use <- NA
                y_use <- NA
            } 
            dat_sub_sub$y_use <- y_use
            dat_sub_sub$window_center <- x_use
            dat_sub_sub[[new_variable_name]] <- dat_sub_sub$rel_log_neg_current - y_use
            dat_sub_sub
        }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .) -> df

    avg_peak_center <- mean(df$window_center, na.rm=TRUE)

    lapply(levels(df$conc_factor), function(level) {
        dat_sub <- subset(df, conc_factor==level)
        lapply(levels(dat_sub$device), function(dev) {
            dat_sub_sub <- subset(dat_sub, device==dev)
            if (is.na(dat_sub_sub$window_center[1])) {
                x_use <- avg_peak_center
                y_use <- approx(dat_sub_sub$potential,
                                dat_sub_sub$rel_log_neg_current,
                                avg_peak_center)$y
                dat_sub_sub$y_use <- y_use
                dat_sub_sub$window_center <- x_use
            } else {
                y_use <- dat_sub_sub$y_use[1]
            }
            dat_sub_sub[[new_variable_name]] <- dat_sub_sub$rel_log_neg_current - y_use
            dat_sub_sub
        }) %>% do.call(rbind, .)
    }) %>% do.call(rbind, .)
}

make_signal_extractor_hessian <- function() {
    function(df) {
        sub_df <- subset(df, in_window)
        sm <- smooth.spline(sub_df$potential, sub_df$resid)
        window_center <- sub_df$window_center[1]
        res_df <- df[1, c("conc_factor", "device", "conc"), drop=FALSE]
        res_df$signal <- as.vector(-Rdistance::secondDeriv(window_center, FUN=function(x) {predict(sm, x)$y}))
        res_df
    }
}

window_data <- function(df, new_variable_name, window_width) {
    lapply(levels(df$conc_factor), function(level) {
        dat_sub <- subset(df, conc_factor==level)
        lapply(levels(dat_sub$device), function(dev) {
            dat_sub_sub <- subset(dat_sub, device==dev)
            window_center <- dat_sub_sub[1, "window_center"]
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
                                         make_residual_calculator,
                                         residual_calculator_params,
                                         make_signal_extractor,
                                         signal_extractor_params,
                                         curr_var,
                                         plot_file_prefix=NULL) {

    ## make the Gaussian Process smoother that works on the windowed votammagram
    residual_calculator <- do.call(make_residual_calculator,
                                   residual_calculator_params)

    ## make the function that extracts the signal
    signal_extractor <- do.call(make_signal_extractor,
                                signal_extractor_params)

    ## make a data frame with the original data, but also with a new column "resid"
    ## containing the residuals from the Gaussian Process fit
    raw_fac_with_resid <- raw_df_fac %>%
        split(list(conc_factor=raw_df_fac$conc_factor, device=raw_df_fac$device)) %>%
        lapply(residual_calculator) 

    ## stack the data frames (one for each voltammogram) into a single tall data frame
    raw_fac_with_resid_df <- do.call(rbind, raw_fac_with_resid)

    ## make a data frame with the conc level, the device number, and the final signal value
    res_with_signals <- raw_fac_with_resid %>%
        lapply(signal_extractor) %>%
        do.call(rbind, .)

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
        plot_conc_faceted_voltammograms(raw_fac_with_resid_df, "resid", "signal",
                                       paste(plot_file_prefix,
                                             "voltammograms-log-rel-sub-resid.pdf",
                                             sep=""))

        p <- ggplot2::ggplot(data=res_with_signals,
                             ggplot2::aes(conc, signal, colour=device)) +
            ggplot2::theme_classic() +
            ggplot2::geom_point()
        ggplot2::ggsave(p, file=paste(plot_file_prefix, "-dot-plot.pdf", sep=""), width=3, height=3)        
    }

    results <- list(avg_rel_err_l1=pct_error_avg_l1,
                    avg_rel_err_l2=pct_error_avg_l2,
                    r=corrcoef,
                    r2=corrcoef*corrcoef,
                    all_results=res_with_signals)

    ## return the performance data via a list
    results
}

#' Run the script `process-voltammograms.R` that is installed in the `ldnh` package.
#' @export
#' @examples
#' process_voltammograms()

process_voltammograms <- function() {
    fork::signal("SIGINT", "default")

    file_stdin <- file('stdin', 'r')

    args <- commandArgs(trailingOnly=TRUE)
    if (length(args) > 0) {
        file_spec <- args[1]
        cat(sprintf("Loading input from file: %s", file_spec))
    } else {

        cat("Specify the file (xlsx format) of melted data to load: ")
        file_spec <- readLines(file_stdin, n=1)

        if (nchar(file_spec)==0) {
            stop("response to prompt was empty; expected a nonzero-length response")
        }
    }

    if (! file.exists(file_spec)) {
        stop(sprintf("File does not exist: %s", file_spec))
    }

    cat("Load worksheet number: [1] ")

    sheet_index <- readLines(file_stdin, n=1)

    if (nchar(sheet_index)==0) {
        sheet_index <- 1
    } else {
        sheet_index <- as.integer(sheet_index)
    }

    cat(sprintf("Specify the starting voltage for the detilting procedure, as a floating-point number [default: %f] ", detilt_start_v_default))
    detilt_start_v <- readLines(file_stdin, 1)
    if (nchar(detilt_start_v) == 0) {
        detilt_start_v <- detilt_start_v_default
    } else {
        detilt_start_v <- as.numeric(detilt_start_v)
    }

    cat(sprintf("Specify the ending voltage for the detilting procedure, as a floating-point number [default: %f] ", detilt_end_v_default))
    detilt_end_v <- readLines(file_stdin, 1)
    if (nchar(detilt_end_v) == 0) {
        detilt_end_v <- detilt_end_v_default
    } else {
        detilt_end_v <- as.numeric(detilt_end_v)
    }


    cat(sprintf("Specify the starting voltage for the peakfinding procedure, as a floating-point number [default: %f] ", peakfind_start_v_default))
    peakfind_start_v <- readLines(file_stdin, 1)
    if (nchar(peakfind_start_v) == 0) {
        peakfind_start_v <- peakfind_start_v_default
    } else {
        peakfind_start_v <- as.numeric(peakfind_start_v)
    }

    cat(sprintf("Specify the ending voltage for the peakfinding procedure, as a floating-point number [default: %f] ", peakfind_end_v_default))
    peakfind_end_v <- readLines(file_stdin, 1)
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
    output_file_spec_prefix <- readLines(file_stdin, 1)
    if (nchar(output_file_spec_prefix) == 0) {
        output_file_spec_prefix <- output_file_spec_prefix_default
    }

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
        
        ## de-tilt the voltammograms by subtracting the average of the blanks
        ##        detilt_data_avg_blanks("rel_log_neg_current") magrittr::%>%
        detilt_linear_model("rel_log_neg_current", c(detilt_start_v, detilt_end_v)) |>

    ## adjust the voltammograms so that they have zero signal value at the "window center";
    ## we are going to search for a local maximum in the range [1.0, 1.1] volts
    make_curves_kiss_avg("rel_log_neg_current_sub",
                               c(peakfind_start_v, peakfind_end_v)) |>
    ## subset the data to a window of 0.1 volts width around the center that we found
    window_data("in_window", 0.5*(peakfind_end_v - peakfind_start_v)) |>
    
    analyze_fitter_and_extractor(
              make_residual_calculator_passthrough,
              residual_calculator_params=list(current_var="rel_log_neg_current_sub"),
              make_signal_extractor=make_signal_extractor_hessian,
              signal_extractor_params=list(),
              curr_var="rel_log_neg_current_sub",
              plot_file_prefix=output_file_spec_prefix) -> results

    xlsx::write.xlsx(results$all_results,
                     file=paste(output_file_spec_prefix, "-processed.xlsx", sep=""),
                     col.names=TRUE,
                     row.names=FALSE)

    xlsx::write.xlsx(data.frame(results[names(results)[which(names(results) != "all_results")]]),
                     file=paste(output_file_spec_prefix, "-summary.xlsx", sep=""),
                     col.names=TRUE,
                     row.names=FALSE)
}
