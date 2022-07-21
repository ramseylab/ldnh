# version 0.1.3.9000

    - Now uses Gaussian kernel smoothing before the detiling procedure,
    for compatibility with data from the SIMStat.
    
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
