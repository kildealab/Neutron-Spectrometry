# Instructions for `unfold_trend.exe`

This application is used to unfold measurements obtained with a Nested Neutron Spectrometer and
record the value of a parameter of interest at each MLEM iteration. Example usage is to output the
average relative error between measured values and MLEM-reconstructed values at each iteration to
a CSV file. These data can then be plotted using `plot_lines.exe`.

## Table of Contents

* [Input files](#input-files)
    * [Measurements file](#measurements-file)
    * [Settings file](#settings-file)
    * [Energy bins](#energy-bins)
    * [NNS response functions](#nns-response-functions)
    * [Guess spectrum](#guess-spectrum)
* [Output files](#output-files)
    * [Trend file](#trend-file)
* [Settings](#settings)

## Input files

### Measurements file
* This file contains the measured NNS data to be unfolded.
* Default file: `input/measurements.txt`
* Can specify an alternative measurements file within the settings file via the `path_measurements` parameter (described below).
* The required format of the file is shown in `input/template_measurements.txt`.
    * Start with 4 header lines to: (1) describe the measurements, (2) specify delivered dose per moderator shell [MU], (3) specify dose rate used [MU/min], (4) specify the irradiation duration [seconds] 
    * Structure of the data portion is such that measured data can be readily copy/pasted from Google Sheets.
* Input values should already have noise and photon contamination removed.
* If multiple measurements were obtained for the same moderator shell, input them sequentially on subsequent lines and adjust the `num_mes_per_shell` setting accordingly.

### Settings file
* This file contains all of the user-configurable settings for the application.
* Default file: `input/unfold_trend.cfg`
* Can specify alternative settings file at runtime via:
```
./unfold_trend.exe --configuration <file_name>
```
* The description of each setting is provided in the [Settings Table below](#settings).
* Default values are indicated where applicable.
    * To use default values, **do not delete settings, simply leave the value blank**.

### Energy bins
* This file contains the energies over which the NNS response functions are defined, and thus the energy bins at which unfolded neutron fluence spectra are quantified.
* Units: [MeV]
* The first value is the lowest energy.
* Values are placed on subsequent lines (no commas).
* File is set via the `path_energy_bins` setting.

### NNS response functions
* This file contains the values of the NNS response functions for all eight moderator configurations over the energy range specified in the [energy bins file](#energy-bins).
* The response of a given moderator is presented on a single line, with values separated by commas and from low to high energy.
* Units: [cm<sup>2</sup>]
* The response of the bare detector is on the first line.
* Subsequent lines have the response for the next moderator configuration.
* File is set via the `path_system_response` setting.

### Guess spectrum
* This file contains a guess neutron fluence spectrum used as the starting point of MLEM unfolding, defined over the energy range specified in the [energy bins file](#energy-bins).
* Units: [neutrons cm<sup>-2</sup> s<sup>-1</sup>]
* The first value corresponds to the lowest energy bin.
* Values are placed on subsequent lines (no commas).
* File is set via the `path_input_spectrum` setting.

## Output files

### Trend file
* This CSV file contains the output of the application.
* The format of the file varies with the `algorithm` used:
    * If `algorithm=mlem`:
        * The first line specifies the iteration indices at which the `parameter_of_interest` was calculated.
        * The second line contains the comma-separated values of the `parameter_of_interest`.
        * If the application is ran multiple times using the same output file, the output data will be appended on new lines to the existing file.
    * If `algorithm=correction_factors`:
        * The first line specifies the energy bins.
        * The first column specifies the iteration indices.
        * The remaining cells comprise a 2D matrix of the spectral correction factors.
    * If `algorithm=trend`:
        * The first line specifies the # of moderators.
        * The first column specifies the iteration indices.
        * The remaining cells comprise a 2D matrix of the ratios between the MLEM reconstructed and measured values.
* File is set via the `path_output_trend` parameter.


## Settings:

| Name | Default value | description |
| ---- | ------------- | ----------- |
| `algorithm` | `mlem` | Specify which unfolding algorithm to use.<br>`mlem`: Calculate `parameter_of_interest` at specified iteration intervals.<br>`map`: Calculate `parameter_of_interest` at specified iteration and beta intervals.<br>`correction_factors`: Calculate correction factors applied to every spectral value at specified iteration intervals.<br>`trend`: Calculate ratio between MLEM-reconstructed values and measurements at specified iteration intervals. |
| `beta_max` | `1E-8` | Applicable if `algorithm=map`. Use in conjunction with `beta_min` to specify range of beta values over which to calculate `parameter_of_interest`. Log-10 intervals are used between min and max values. |
| `beta_min` | `1E-10` | Applicable if `algorithm=map`. Use in conjunction with `beta_max` to specify range of beta values over which to calculate `parameter_of_interest`. Log-10 intervals are used between min and max values. |
| `derivatives` | `0` | When using `algorithm=mlem`, set `derivatives=1` if want to calculate the rate of change of change (derivative) of the `parameter_of_interest`. |
| `f_factor` | `7.2` | Conversion coefficient between neutron current and CPS for NNS [fA/cps]. |
| `iteration_increment` | `100` | Use in conjuction with `iteration_min` and `iteration_max` to specify the range and intervals of MLEM iterations at which to calculate parameter of interest. E.g. start at 100 iterations, calculate every 100 iterations, and stop at 10,000 iterations. Value must evenly divide the difference between min and max. |
| `iteration_max` | `10000` | See `iteration_increment`. |
| `iteration_min` | `100` | See `iteration_increment`. |
| `meas_units` | `nc` |  Specify units of measured values {`nc`,`cps`}. |
| `nns_normalization` | `1.14` | NNS-dependent normalization factor. |
| `num_meas_per_shell` | `1` | # of measured values input per moderator shell. |
| `parameter_of_interest` | `total_fluence` | Parameter to be calculated at specified iterations. {`avg_mlem_ratio`,`chi_squared_g`,`j_factor`,`j_factor2`,`max_mlem_ratio`,`noise`,`nrmsd`,`reduced_chi_squared`,`rms`,`total_dose`,`total_fluence`} |
| `path_energy_bins` | `input/energy_bins.csv` | Pathname to [energy bins file](#energy-bins). |
| `path_icrp_factors` | `input/`<br>`icrp_conversion_coefficients` | Pathname to [file containing ambient dose equivalent conversion coefficients](#ambient-dose-equivalent-conversion-factors) [pSv cm^2]. |
| `path_input_spectrum` | `input/spectrum_step.csv` | Pathname to [input (guess) spectrum file](#guess-spectrum). |
| `path_measurements` | `input/measurements.txt` | Pathname to [NNS measurements file](#measurements-file). |
| `path_output_trend` | `output/output_trend.csv` | Pathname to CSV file to store output trend. |
| `path_ref_spectrum` | N/A | Pathname to a spectrum file to be used as ground-truth reference spectrum when `parameter_of_interest=rms`. |
| `path_system_response` | `input/response_nns_he3.csv` | Pathname to [NNS response functions file](#nns-response-functions). |
| `prior` | `mrp` | Type of prior calculation to be done if `algorithm=map`.<br>`quadratic`: smoothing, no edge preservation.<br>`mrp`: median root prior; preserves edges by not penalizing regions of monotonic increase or decrease.<br>`medianrp`: mean root prior; custom written; similar to `mrp` but based on mean of neighbours. |
| `trend_type` | `cps` | Use if `algorithm=trend`. Defines how first output row containing measured values appears.<br>`ratio`: all will be 1 (ratio with itself).<br>`cps`: output the measured values in CPS. |