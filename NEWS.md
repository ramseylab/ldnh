# version 0.1.9.9000

    - Adding new column `avg_cv` to the `-summary.xlsx` spreadsheet

# version 0.1.8.9000

    - Fixing problem with `convert_raw_to_melted_xlsx` not working when the
    raw potentiostat data file has extra lines in its file header
    
# version 0.1.7.9000

    - Changing default voltage for the right side of the normalization (i.e.,
    "detilting") window to 1.135.
    - Improving plotting code.
   
# version 0.1.6.9000

    - Fixing bug with export of `process_voltammograms`

# version 0.1.5.9000

    - Fixing bug where the `stringr` dependency wasn't listed in DESCRIPTION.
    
# version 0.1.4.9000

    - Adding Gaussian kernel smoothing between the log-transformation step
    and the detilting step, for compatibility with voltammograms from the
    SIMStat.
    
# version 0.1.3.9000

    - The `process_voltammograms` function now will not give an error if there
    are extra named (or unnamed/empty) columns in the input spreadsheet of 
    melted data, other than the four required columns named `potential`,
    `device`, `conc`, and `current`. The script will now just ignore the
    extra columns.

# version 0.1.2.9000

    - Script can now handle the case where there are different combinations
    of device identifiers (i.e., device/replicate numbers) used at different
    analyte concentration levels. Prior to this fix, in such a use-case, 
    the script would error out with an obscure and un-helpful error message.

# version 0.1.1.9000

    - In the `-processed.xlsx` file, save (for each sample) the potential at
      which the local maximum (analyte peak) was found, in column
      `peak_potential`.
    
# version 0.1.0.9000

    - Change analysis method to use a smoothed spline fit to censored data, for
      background model
    - Added MIT license.
    
# version 0.0.2.9000

    - Merged to `development` branch on 2022-07-14
	- Provides the `convert_raw_to_melted_xlsx` function

# version 0.0.1.9000

    - Merged to `main` branch on 2022-07-14
	- First version released; provides the `process_voltammograms` function
