To install this R package, run Rstudio and then follow these steps, in order:

1. Install the R package `remotes` from CRAN by running the following command in the Rstudio console:

```install.packages("remotes")```

2. Install the R package `ldnh` from GitHub by running the following command in the Rstudio Console:

```remotes::install_github("ramseylab/ldnh")```

You may be asked "These packages have more recent versions
available; which would you like to update?"; in that case,
respond `3` (corresponding to "None").

You may be asked if you wish to install packages from source; respond "no".

3. In the Rstudio Console, run the following command in order to run the analysis script:

```ldnh::process_voltammograms()```

Some notes:

- Once you have performed Step 1 once, you should not need to perform it again.
- Once you have performed Step 2 once, you should only need to perform it again
unless there is an update to the `ldnh` R package that you want to use within
your local Rstudio installation.
- You need to perform Step 3 above, each time you want to process some voltammograms.

