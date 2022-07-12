
# Requirements

You will need R in order to use this software. You can get R for free at the 
[CRAN website](https://cran.r-project.org/).

# Installing Rstudio

If you do not have Rstudio installed, you might want to install it; Rstudio 
provides an Integrated Development Environment (IDE) for using R. You can get
Rstudio for free at the
[Rstudio website](https://www.rstudio.com/products/rstudio/download/).

Rstudio is not required in order to use the `ldnh` R package.

# Installing the `ldnh` package into R

To install this R package, run Rstudio (or R) and then follow these steps, in order:

1. Install the R package `remotes` from CRAN by running the following command in the Rstudio console:

```install.packages("remotes")```

2. Install the R package `ldnh` from GitHub by running the following command in the R Console:

```
options(install.packages.compile.from.source = "never")
remotes::install_github("ramseylab/ldnh")
```

You may be asked "These packages have more recent versions
available; which would you like to update?"; in that case,
respond `3` (corresponding to "None").

3. If you have updated the `ldnh` package (versus installing it for the first time),
you should quit and restart R; otherwise, just proceed to the section "Running the
`ldnh` software to process voltammograms" below.

# Running the `ldnh` software to process voltammograms

In the R console, run the following command in order to run the analysis script:

```ldnh::process_voltammograms()```

Some notes:

- Once you have performed Step 1 once, you should not need to perform it again.
- Once you have performed Step 2 once, you should only need to perform it again
unless there is an update to the `ldnh` R package that you want to use within
your local R installation.
- You need to perform Step 3 above, each time you want to process some voltammograms.

# Input spreadsheet format:

The script `process_voltammograms` (see Step 3 above) expects to load your voltammogram
data in an Excel spreadsheet (an `.xlsx` file) with the following four columns in it.
The columns do not have to be in any particular order, but their names must be 
*exactly* as shown below:

- `potential`: this column contains the voltage levels (numeric type)
- `current`: this column contains the peak current levels, at each potential level (numeric type)
- `device`: identifies the specific device (i.e., replicate) (integer type)
- `conc`: identifies the concentration of analyte (numeric type)

As you can see from the above, the spreadsheet is required to have the data
arranged in a "melted" format of four-tuples, like this:

| potential | device | conc | current |
| --------- | ------ | ---- | ------- |
| 0.504     | 1      | 0    | \-4.858 |
| 0.508     | 1      | 0    | \-4.758 |
| 0.512     | 1      | 0    | \-4.785 |
| 0.516     | 1      | 0    | \-4.825 |
| 0.52      | 1      | 0    | \-4.88  |
| 0.524     | 1      | 0    | \-4.931 |
| 0.528     | 1      | 0    | \-4.99  |
| 0.532     | 1      | 0    | \-5.04  |
| 0.536     | 1      | 0    | \-5.097 |
| 0.54      | 1      | 0    | \-5.151 |
| 0.544     | 1      | 0    | \-5.206 |
| 0.548     | 1      | 0    | \-5.26  |
| 0.552     | 1      | 0    | \-5.314 |
| 0.556     | 1      | 0    | \-5.367 |
| 0.56      | 1      | 0    | \-5.424 |
| 0.564     | 1      | 0    | \-5.477 |
| 0.568     | 1      | 0    | \-5.539 |
| 0.572     | 1      | 0    | \-5.589 |
| 0.576     | 1      | 0    | \-5.645 |
| 0.58      | 1      | 0    | \-5.707 |
| 0.584     | 1      | 0    | \-5.751 |
| 0.588     | 1      | 0    | \-5.81  |
| 0.592     | 1      | 0    | \-5.873 |
| 0.596     | 1      | 0    | \-5.927 |
| 0.6       | 1      | 0    | \-5.983 |
| 0.604     | 1      | 0    | \-6.044 |
| 0.608     | 1      | 0    | \-6.101 |
| 0.612     | 1      | 0    | \-6.162 |
| 0.616     | 1      | 0    | \-6.228 |
| 0.62      | 1      | 0    | \-6.282 |
| 0.624     | 1      | 0    | \-6.343 |
| 0.628     | 1      | 0    | \-6.407 |
| 0.632     | 1      | 0    | \-6.472 |
| 0.636     | 1      | 0    | \-6.536 |
| 0.64      | 1      | 0    | \-6.598 |
| 0.644     | 1      | 0    | \-6.662 |
| 0.648     | 1      | 0    | \-6.726 |
| 0.652     | 1      | 0    | \-6.794 |

