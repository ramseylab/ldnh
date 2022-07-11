#'
#' @import xlsx remotes Rdistance magrittr ggplot2 fork SplinesUtils
#'
function() {}  ## this is required in order to force roxygen2::roxygenise to
               ## generate "import" directives in the NAMESPACE file

#' Run the script `process-voltammograms.R` that is installed in the `ldnh` package.
#' @export
#' @examples
#' process_voltammograms()

process_voltammograms <- function() {
    source(paste(.libPaths(), "ldnh/exec/process-voltammograms.R", sep="/"))
}
