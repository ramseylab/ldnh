To install this R package, run Rstudio and then follow these steps, in order:

1. Install the R package `remotes` from CRAN:

```install.package("remotes")```

2. Install the R package `ldnh` from GitHub:

```remotes::install_github("ramseylab/ldnh")```

You may be asked "These packages have more recent versions
available; which would you like to update?"; in that case,
respond `3` (corresponding to "None").

You may be asked if you wish to install packages from source; respond "no".

3. Run the following code, in order to run the analysis script:

```source(paste(.libPaths(), "ldnh/exec/process-voltammograms.R", sep="/"))```


