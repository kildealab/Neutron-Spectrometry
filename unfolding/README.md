# **Instructions** #
README.md

Authors: Logan Montgomery, Georges Al Makdessi, Robert Maglieri, John Kildea

This file explains:
1) General instructions for compiling and executing a program
2) The purpose of the programs
3) The file structure of the project, including inputs & outputs

Please refer to dependencies.txt for dependencies

### Execution Instructions ###
* Compile the program:
    $ make

* Execute your program of choice. For example, if unfolding, run the following
    $ ./unfold_spectrum.exe
    * Additional options may be provided at the command line for unfolding. These options and the default values are:
        * --configuration [mlem.cfg]

* ***NOTE*** Various input files, including that containing measured values, must be specified in the configuration file

### Program #1 unfold_spectrum.exe ###
* Unfolds measurements obtained with a Nested Neutron Spectrometer into a neutron fluence spectrum.
* ***Inputs***:
    * Configuration settings (e.g. algorithm, parameter of interest, max number of iterations, minimum number of iterations, step size, etc.)
        * --configuration [auto.cfg]
* ***Outputs***:
    * Unfolded fluence spectrum & Poisson uncertainty are appended as CSV data to a user-specified output file.
* ***Details***:
    * The program can unfold using a standard MLEM algorithm, or using the MAP (or OSL) algorithm with a user-specified prior (specify in map.cfg). Allowed priors are:
        * ***quadratic*** - First prior for MAP/OSL unfolding, proposed by Green, 1990. Smooths data, no edge preservation.
        * ***mrp*** - Median Root Prior used to penalize areas that are not monotonically increasing or decreasing (preserves edges).
* ***Allowed Settings***:
    * mlem_cutoff= Maximum number of iterations
    * nns_normalization= Normalization factor for the specific NNS used for measurements
    * mlem_max_error= Threshold error in ratio between measurements and reconstructed measurements, below which MLEM is stopped
    * f_factor= Conversion coefficient between neutron current and CPS [fA/cps]
    * uncertainty_type= Specify the type of spectral uncertainty calculation {poisson,j_bounds}
    * num_uncertainty_samples= Number of samples to be taken to generate uncertainty region (if poisson/gaussian)
    * num_meas_per_shell= Number of measured values input per moderator shell
    * sigma_j= The uncertainty in J to be used to generate uncertainty region (if j_bounds)
    * meas_units= Specify the units used in the measurements files {nc,cps}
    * algorithm= Specify which algorithm to use {e.g. mlem, map, or mlemstop}
    * beta= Beta factor to use if unfolding using the MAP-EM algorithm
    * prior= Prior to use if unfolding using the MAP-EM algorithm {mrp,quadratic}
    * cps_crossover= Crossover CPS value to use in calculating J_threshold for MLEM-STOP
    * generate_report= 1 or 0 to indicate whether a report should be generated for the unfolding process
    * generate_figure= 1 or 0 to indicate whether the unfolded spectrum should be plotted
    * path_output_spectra = Path and filename (relative or absolute) of the file for unfolded spectra
    * path_report = Path and filename (relative or absolute) of the file for the unfolding report
    * path_figure = Path and filename (relative or absolute) of the file for the plotted spectrum
    * measurements_path= Path and filename (relative or absolute) of the measurments file
    * input_spectrum_path= Path and filename (relative or absolute) of the input/guess spectrum
    * energy_bins_path= Path and filename (relative or absolute) of the energy bins of the detector
    * system_response_path= Path and filename (relative or absolute) of the detector response functions
    * icrp_factors_path= Path and filename (relative or absolute) of the ICRP ambient dose equivalent conversion coefficients

### Program #2 plot_spectra.exe ###
* Plots one or more neutron spectra (and their uncertainties) on a single set of axes.
* ***Inputs***:
    * Plot settings are read in from input/plot.cfg. Please see below for an explanation of the settings.
        * Note, settings may be left empty and default values will be assigned
    * Spectra are read in from the CSV input file specified in input/plot_spectra.cfg (default is output/output_spectra.csv)
* ***Outputs***:
    * Spectra are plotted onto a file specified in input/plot.cfg (default is output/output_spectra.png)
* ***Details***:
    * The program is designed to plot neutron spectra that may be easily modified using the settings specified by the user in input/plot.cfg. The program must be rerun when the settings are modified, but does not require recompilation!
* ***Allowed Settings***:
    * input_filename= filename including extension
    * input_dir= pathname ending with / (e.g. input/)
    * output_filename= filename including extension
    * output_dir= pathname ending with / (e.g. output/)

    * title= Title displayed at the top of the 
    * x_label= x-axis title
    * y_label= y-axis title

    * x_res= horizontal resolution (dimension) as a number (e.g. 3200)
    * y_res= vertical resolution (dimension) as a number (e.g. 2400)

    * x_min= minimum x-axis value
    * x_max= maximum x-axis value
    * y_min= minimum y-axis value
    * y_max= maximum y-axis value
    * y_num_divs= Number that specifies number of major and minor divisions (tick marks). See https://root.cern.ch/doc/master/classTAttAxis.html#ae3067b6d4218970d09418291cbd84084
    * y_digits_max = Max number of digits allowed on the y axis labels before scientific notation is used

    * color_series= Comma-delimited list of HEX color values (e.g. #000000,#FFFFFF)
    * grayscale= 1 or 0 indicating whether figure should be output in grayscale
    * line_style= Comma-delimited list of numeric line styles. See https://root.cern.ch/doc/master/classTAttLine.html
    * line_width= Commad-delimited list of numeric line widths. See https://root.cern.ch/doc/master/classTAttLine.html
    * border_width= Numeric width for the border (axes) of the graph

    * show_error= Comma-delimited list of 1 or 0 indicating whether uncertainty should be displayed for each spectrum
    * error_style= The visual style with which uncertainties are displayed. See https://root.cern.ch/doc/master/classTHistPainter.html#HP01b
    * error_fill_style= The fill pattern to use in the uncertainty region. See https://root.cern.ch/doc/master/classTAttFill.html
    * rows_per_spectrum= The number of rows in the input file dedicated to each spectrum (default is 3: 1 for the spectrum, 1 for upper uncertainty, and 1 for lower uncertainty)

    * legend= 1 or 0 indicating whether to include the legend
    * legend_entries= Comma-delimited list of text entries to put in the legend
    * legend_coords= 4 comma-delimited double values representing coordinates of the legend. Startx, Starty, Endx, Endy.

    * textbox= 1 or 0 indicating whether to include a textbox
    * textbox_coords= 4 comma-delimited double values representing coordinates of the legend. Startx, Starty, Endx, Endy.
    * textbox_text= Comma-delimited list of text entries to put in the textbox

    * plot_per_mu= 1 or 0 indicating whether to plot per MU instead of per second
    * number_mu= Comma-delimited list of number of MU delivered for each spectrum; required to plot per MU
    * duration= Comma-delimited list of duration of delivery for each spectrum; required to plot per MU

    * normalize= 1 or 1 indicating whether each spectrum should be normalized such that its max value is set to 1

### Program #3 auto_unfold_spectrum.exe ###
* Unfold a set of neutron measurements and output a CSV file containing values of a parameter of interest that was calculated iteratively during the unfolding process.
* ***Inputs***:
    * Configuration settings (e.g. algorithm, parameter of interest, max number of iterations, minimum number of iterations, step size, etc.)
        * --configuration [auto.cfg]
* ***Outputs***:
    * A CSV file containting the parameter of interest values are output to a file specified in input/auto.cfg.
        * For MLEM, POI values are calculated as a function of N (one dimension), and are appended to the output file if it previously contained data.
        * For MAP, POI values are calculated as a function of N and beta (two dimensions), and overwrite the output file if it previously contained data.
* ***Details***:
    * The output files are intended to be readily used with the plot_lines.exe program (for MLEM data) and the plot_surface.exe program (for MAP data).
    * Allowed parameter of interest:
        * total_fluence
        * total_dose
        * total_energy_correction (MAP specific)
        * max_mlem_ratio
        * avg_mlem_ratio
* ***Allowed Settings***:
    * nns_normalization= Normalization factor for the specific NNS used for measurements
    * f_factor= Conversion coefficient between neutron current and CPS [fA/cps]
    * meas_units= Specify the units used in the measurements files {nc,cps}
    * beta= Beta factor to use if unfolding using the MAP-EM algorithm
    * prior= Prior to use if unfolding using the MAP-EM algorithm {mrp,quadratic}
    * min_num_iterations= First iteration at which to score the parameter of interest
    * max_num_iterations= Last iteration at which to score the parameter of interest
    * iteration_increment= How many iterations to increment between scoring parameter of interest
    * min_beta= Minimum beta at which to score parameter of interest for MAP-EM
    * max_beta= Maximum beta at which to score parameter of interest for MAP-EM
    * parameter_of_interest= The parameter of interest to be scored if using the mlem algorithm {rms,total_fluence,total_dose,max_mlem_ratio,avg_mlem_ratio,j_factor}
    * derivatives= 1 of 0 indicating if you want to score the derivative of the parameter of interest with respect to the iteration index
    * auto_output_path= Path and filename (relative or absolute) of the output (results) file
    * algorithm= Specify what algorithm to use {mlem,map,trend,correction_factors}
    * trend_type= Use if alogirhtm = 
    * measurements_path= Path and filename (relative or absolute) of the measurments file
    * input_spectrum_path= Path and filename (relative or absolute) of the input/guess spectrum
    * energy_bins_path= Path and filename (relative or absolute) of the energy bins of the detector
    * system_response_path= Path and filename (relative or absolute) of the detector response functions
    * icrp_factors_path= Path and filename (relative or absolute) of the ICRP ambient dose equivalent conversion coefficients
    * ref_spectrum_path= Path and filename (relative or absolute) of the reference spectrum used with RMSE calculations


### Program #4 plot_lines.exe ###
* Plots one or more data series (y as a function of x) on a single set of axes.
* ***Inputs***:
    * Plot settings are read in from input/plot_lines.cfg. Please see below for an explanation of the settings.
        * Note, settings may be left empty and default values will be assigned
    * Data are read in from the CSV input file specified in input/plot_lines.cfg (default is output/poi_output_mlem.csv)
* ***Outputs***:
    * Plots are created in a file specified in input/plot_lines.cfg (default is output/poi_output_mlem.png)
* ***Details***:
    * The program is designed to generate plots of parameter of interest values that may be easily modified using the settings specified by the user in input/plot_lines.cfg. The program must be rerun when the settings are modified, but does not require recompilation!
* ***Allowed Settings***:
    * input_filename= filename including extension
    * input_dir= pathname ending with / (e.g. input/)
    * output_filename= filename including extension
    * output_dir= pathname ending with / (e.g. output/)

    * data_format= Specify the format of the data rows read in from CSV to be plotted {xyy,xyxy}

    * title= Title displayed at the top of the 
    * x_label= x-axis title
    * y_label= y-axis title
    * x_label_offset= Double value to shift the x-axis title up or down (1.0 is the default)
    * y_label_offset= Double value to shift the y-axis title up or down (1.0 is the default)

    * x_res= horizontal resolution (dimension) as a number (e.g. 3200)
    * y_res= vertical resolution (dimension) as a number (e.g. 2400)

    * x_min= minimum x-axis value
    * x_max= maximum x-axis value
    * y_min= minimum y-axis value
    * y_max= maximum y-axis value
    * x_num_divs= Number that specifies number of major and minor divisions (tick marks). See https://root.cern.ch/doc/master/classTAttAxis.html#ae3067b6d4218970d09418291cbd84084
    * y_num_divs= Number that specifies number of major and minor divisions (tick marks). See https://root.cern.ch/doc/master/classTAttAxis.html#ae3067b6d4218970d09418291cbd84084

    * x_log= 1 or 0 indicating whether the x axis should be logarithmic
    * y_log= 1 or 0 indicating whether the y axis should be logarithmic

    * x_grid= 1 or 0 indicating whether there should be grid lines oriented parallel to the x axis
    * y_grid= 1 or 0 indicating whether there should be grid lines oriented parallel to the y axis

    * plot_type= Character representing the plot type. Allowed options are l (for line plots) and p (for scatter plots)
    * color_series= Comma-delimited list of HEX color values (e.g. #000000,#FFFFFF)
    * grayscale= 1 or 0 indicating whether figure should be output in grayscale
    * line_style= Comma-delimited list of numeric line styles. See https://root.cern.ch/doc/master/classTAttLine.html
    * line_width= Commad-delimited list of numeric line widths. See https://root.cern.ch/doc/master/classTAttLine.html
    * border_width= Numeric width for the border (axes) of the graph
    * marker_style= Comma-delimited list of numeric marker styles. See https://root.cern.ch/doc/master/classTAttMarker.html
    * marker_size= Comma-delimited list numeric marker sizes. See https://root.cern.ch/doc/master/classTAttMarker.html

    * margin_left= Double value representing the fraction of the plot size composed of a left margin
    * margin_right= Double value representing the fraction of the plot size composed of a right margin
    * margin_top= Double value representing the fraction of the plot size composed of a top margin
    * margin_bottom= Double value representing the fraction of the plot size composed of a bottom margin

    * legend= 1 or 0 indicating whether to include a legend
    * legend_entries= Comma-delimited list of text entries to put in the legend
    * legend_coords= 4 comma-delimited double values representing coordinates of the legend. Startx, Starty, Endx, Endy.
    * legend_border_size= Integer representing legend border size
    * legend_text_size= Double value representing text size as fraction of plot size (default is 0.035)

    * textbox= 1 or 0 indicating whether to include a textbox
    * textbox_coords= 4 comma-delimited double values representing coordinates of the legend. Startx, Starty, Endx, Endy.
    * textbox_text= Comma-delimited list of text entries to put in the textbox

### Program #5 plot_surface.exe ###
* Plots z values as a function of x and y.
* ***Inputs***:
    * Plot settings are read in from input/plot_surface.cfg. Please see below for an explanation of the settings.
        * Note, settings may be left empty and default values will be assigned
    * Data are read in from the CSV input file specified in input/plot_surface.cfg (default is output/poi_output_mlem.csv)
* ***Outputs***:
    * Plots are created in a file specified in input/plot_lines.cfg (default is output/poi_output_mlem.png)
* ***Details***:
    * The program is designed to generate plots of parameter of interest values that may be easily modified using the settings specified by the user in input/plot_lines.cfg. The program must be rerun when the settings are modified, but does not require recompilation!
* ***Allowed Settings***:
    * input_filename= filename including extension
    * input_dir= pathname ending with / (e.g. input/)
    * output_filename= filename including extension
    * output_dir= pathname ending with / (e.g. output/)

    * title= Title displayed at the top of the 
    * x_label= x-axis title
    * y_label= y-axis title
    * z_label= z-axis title

    * x_res= horizontal resolution (dimension) as a number (e.g. 3200)
    * y_res= vertical resolution (dimension) as a number (e.g. 2400)

    * x_min= minimum x-axis value
    * x_max= maximum x-axis value
    * y_min= minimum y-axis value
    * y_max= maximum y-axis value
    * z_min= minimum z-axis value
    * z_max= maximum z-axis value
    * x_num_divs= Number that specifies number of major and minor divisions (tick marks). See https://root.cern.ch/doc/master/classTAttAxis.html#ae3067b6d4218970d09418291cbd84084
    * z_num_divs= Number that specifies number of major and minor divisions (tick marks). See https://root.cern.ch/doc/master/classTAttAxis.html#ae3067b6d4218970d09418291cbd84084

    * color_palette= Number that specifies the ROOT color palette to use. See "High quality predefined palettes" at https://root.cern.ch/doc/master/classTColor.html
    * num_color_bins= The number of different color bins. More bins = more finely resolved surface.
    * border_width= Numeric width for the border (axes) of the graph

### Notes on input files###
* Files read by programs are included in input/
    * ***Energy bins*** - Energy bins corresponding to the values in other input files 
        * Units: MeV
        * Default: energy_bins.csv 

    * ***Input spectrum*** - The "guess" input neutron flux spectrum inserted into the unfolding algorithm. Currently use a step function (high value at low energies and low value at high energies, step at thermal)
        * Units:[n cm^-2 s^-1]
        * Default: spectrum_step.csv

    * ***NNS Response function*** - Vendor-provided NNS relative response as a function of energy for each level of moderation. Line 1: response for bare probe, Line 2: response for 1 moderator, etc. 
        * Units [cm^2]
        * Default: response_nns_he3.csv

    * ***ICRP conversion factors*** - Factors to convert neutron flux [n cm^-2 s^-1] to neutron ambient dose equivalent rate [pSv s^-1]. Values were linearly interpolated from tabulated values in ICRP 74 Table A.42 (1st data column).
        * Units: pSv cm^2
        * Default: icrp_conversions.crp