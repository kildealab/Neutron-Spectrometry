# **Instructions for Unfolding** #
README.md
Creation Date:  2017-08-21
Last Modified:  2017-11-20
Authors: Georges Al Makdessi, John Kildea, Robert Maglieri, Logan Montgomery

This file explains:
1) Instructions for executing the program
2) The file structure of the project, including inputs & outputs

Please refer to dependencies.txt for dependencies


### 1) Execution instructions ###
* Input measured data:
    * Create input file (e.g. input/measurements.csv)
    * Refer to input/template_measurements.csv for structure of the file
        * Line 1: String representing the irradiation conditions
        * Line 2: Measurement duration in seconds
        * Lines 3 through 10: Measured charge values in nC (decreasing number of moderators)
            * Negative charge values are acceptable

* Determine which algorithm to use for unfolding:
    * Settings for each allowed algorithm are provided in the input/ directory as \*.cfg files

* Adjust settings (optional):
    * edit the appropriate input/\*.cfg file
    * If necessary, create new input input/\*.csv files (e.g. a different input spectrum)

* Compile the program:
    $ make

* Execute the program and specify 1 or more input files using the options below (if any option is not
  included then a default file is used, as specified in square brackets below):
    $ ./unfold_spectrum.exe
    * Allowed Options:
        * --configuration [mlem.cfg]
        * --measurements [measurements.txt]
        * --energy-bins [energy_bins.csv]
        * --input-spectrum [spectrum_step.csv]
        * --icrp-factors [icrp_conversions.csv]
        * --nns-response [nns_response.csv]


### 2) File structure ###
* Main unfolding code: unfold_spectrum.cpp

* Settings files: input/\*.cfg, 
    * ***mlem_cutoff*** - max # of iterations of unfolding algorithm
    * ***nns_normalization*** - Vendor-provided normalization factor for the NNS used (varies with NNS)
    * ***mlem_max_error*** - The maximum allowed error in reconstructed measurement values from the true measurment values (in nC). Calculated for each unfolding iteration for each measured value (8). Once all values are within this max error, unfolding algorithm will stop.
    * ***f_factor*** - Vendor-provided calibration factor to convert measured current (i.e. measured charge divided by duration) to counts per second (cps). Value in units of [fA / cps]. Value was verified experimentally.
    * ***num_poisson_samples*** - Number of times to repeat unfolding using poisson-sampled data to generate statistical uncertainties.
    * ***beta*** - [MAP only] Tuning (weighting) parameter to set the strength of the MAP/OSL prior.
    * ***prior*** - [MAP only] The prior to be used as the penalty factor in MAP unfolding. Possible values:
        * ***quadratic*** - First prior for MAP/OSL unfolding, proposed by Green, 1990. Smooths data, no edge preservation.
        * ***mrp*** - Median Root Prior used to penalize areas that are not monotonically increasing or decreasing (preserves edges).
    * ***parameter_of_interest*** - [OPTIMIZE only] The parameter of interest to be calculated in order to determine the optimal beta. Possible values:
        * ***total_fluence*** - The integrated fluence of the unfolded spectrum.
        * ***total_dose*** - The total neutron ambident dose equivalent.
        * ***total_energy_correction*** - The total enery correction term (prior \* beta), summed across all bins.
        * ***max_mlem_ratio*** - The maximum deviation in the MLEM ratio of unfolded data relative to measured data.
    * ***min_num_iterations*** - [OPTIMIZE only] Minimum number of iterations to consider
    * ***max_num_iterations*** - [OPTIMIZE only] Maximum number of iterations to consider
    * ***iteration_increment*** - [OPTIMIZE only] When performing optimization, how many iterations to do before recalculating parameter of interest.
    * ***min_beta*** - [OPTIMIZE only] Minimum beta factor to consider
    * ***max_beta*** - [OPTIMIZE only] Maximum beta factor to consider

* Files read by the program are included in input/
    * 4 files from this directory are read by the program, but multiple options may exist
        * ***Energy bins*** - Energy bins corresponding to the values in other input files 
            * Units: MeV
            * Default: energy_bins.csv 

        * ***Input spectrum*** - The "guess" input neutron flux spectrum inserted into the unfolding algorithm. Currently use a step function (high value at low energies and low value at high energies, step at thermal)
            * Units:[n cm^-2 s^-1]
            * Default: spectrum_step.csv

        * ***NNS Response function*** - Vendor-provided NNS relative response as a function of energy for each level of moderation. Line 1: response for bare probe, Line 2: response for 1 moderator, etc. 
            * Units [cm^2]
            * Default: nns_response.csv

        * ***ICRP conversion factors*** - Factors to convert neutron flux [n cm^-2 s^-1] to neutron ambient dose equivalent rate [pSv s^-1]. Values were linearly interpolated from tabulated values in ICRP 74 Table A.42 (1st data column).
            * Units: pSv cm^2
            * Default: icrp_conversions.crp

* Files written to by the program are located in output/
    * ***output_spectrum.csv*** - The unfolded spectrum generated from executing this program is appended to this file. Energy bins are also included. The new spectrum is appended as a column (i.e. a new element is added to each row). Format this way for direct import into spreadsheet program
    * ***output_dose.csv*** - The measured total neutron ambient dose equivalent rate and its uncertainty are appended to this file. The new values are inserted as a new row.
    * ***report_\*.txt*** - Report file generated each time the program is executed. The report contains all pertinent inputs & outputs of the program.
    * ***figure_\*.png*** - Figure file generated each time the program is executed. The figure displays the spectrum and its uncertainty plotted as a function of energy.
    * ***map_file.csv*** - [OPTIMIZE only] 2D "matrix" of parameter of interest values to be used to generate a contour map.

