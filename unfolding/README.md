# **Instructions** #
README.md

Authors: Georges Al Makdessi, John Kildea, Robert Maglieri, Logan Montgomery

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
    * Several additional options may be provided at the command line for unfolding. These options and the default values are:
        * --configuration [mlem.cfg]
        * --measurements [measurements.txt]
        * --energy-bins [energy_bins.csv]
        * --input-spectrum [spectrum_step.csv]
        * --icrp-factors [icrp_conversions.csv]
        * --nns-response [nns_response.csv]

* ***NOTE*** You must specify the measurement data (in measurements.txt or equivalent) and necessary settings (in mlem.cfg or equivalent) before unfolding the data.

### Program #1 unfold_spectrum.exe ###
* Unfolds measurements obtained with a Nested Neutron Spectrometer into a neutron fluence spectrum.
* ***Inputs***:
    * The measured charge data in nC aquired with NNS (as well as dose, dose rate, and measurement duration)
        * --measurements [measurements.txt]
    * Configuration settings for unfolding (e.g. algorithm, stopping criterion, etc.)
        * --configuration [mlem.cfg]
    * Energy bins at which fluence is calculated
        * --energy-bins [energy_bins.csv]
    * Input guess fluence spectrum
        * --input-spectrum [spectrum_step.csv]
    * NNS response functions
        * --nns-response [nns_response.csv]
    * ICRP 74 neutron fluence to neutron ambient dose equivalent conversion coefficients
        * --icrp-factors [icrp_conversions.csv]
* ***Outputs***:
    * Unfolded fluence spectrum & uncertainty are appended as CSV data to output/output_spectra.csv
    * The calculated ambient dose equivalent & uncertainty are appended to output/output_dose.csv
* ***Details***:
    * The program can unfold using a vanilla MLEM algorithm, or using the MAP (or OSL) algorithm with a user-specified prior (specify in map.cfg). Allowed priors are:
        * ***quadratic*** - First prior for MAP/OSL unfolding, proposed by Green, 1990. Smooths data, no edge preservation.
        * ***mrp*** - Median Root Prior used to penalize areas that are not monotonically increasing or decreasing (preserves edges).

### Program #2 plot_spectra.exe ###
* Plots one or more neutron spectra (and their uncertainties) on a single set of axes.
* ***Inputs***:
    * Plot settings are read in from input/plot.cfg. Please refer to input/template_plot.cfg for full list of customizable settings.
    * Spectra are read in from the CSV input file specified in input/plot.cfg (default is output/output_spectra.csv)
* ***Outputs***:
    * Spectra are plotted onto a file specified in input/plot.cfg (default is output/output_spectra.png)
* ***Details***:
    * The program is designed to plot neutron spectra that may be easily modified using the settings specified by the user in input/plot.cfg. The program must be rerun when the settings are modified, but does not require recompilation!

### Program #3 auto_unfold_spectrum.exe ###
* Unfold a set of neutron measurements and output a CSV file containing values of a parameter of interest that was calculated iteratively during the unfolding process.
* ***Inputs***:
    * The measured charge data aquired with NNS (as well as dose, dose rate, and duration)
        * --measurements [measurements.txt]
    * Configuration settings (e.g. algorithm, parameter of interest, max number of iterations, minimum number of iterations, step size, etc.)
        * --configuration [auto.cfg]
    * Energy bins at which fluence is calculated
        * --energy-bins [energy_bins.csv]
    * Input guess fluence spectrum
        * --input-spectrum [spectrum_step.csv]
    * NNS response functions
        * --nns-response [nns_response.csv]
    * ICRP 74 neutron fluence to neutron ambient dose equivalent conversion coefficients
        * --icrp-factors [icrp_conversions.csv]
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

### Program #4 plot_lines.exe ###
* Plots one or more data series (y as a function of x) on a single set of axes.
* ***Inputs***:
    * Plot settings are read in from input/plot_lines.cfg. Please see below for full list of settings.
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

    * plot_type= Character representing the plot type. Allowed options are l (for line plots) and p (for scatter plots)
    * color_series= Comma-delimited list of HEX color values (e.g. #000000,#FFFFFF)
    * line_style= Comma-delimited list of numeric line styles. See https://root.cern.ch/doc/master/classTAttLine.html
    * line_width= Commad-delimited list of numeric line widths. See https://root.cern.ch/doc/master/classTAttLine.html
    * marker_style= Comma-delimited list of numeric marker styles. See https://root.cern.ch/doc/master/classTAttMarker.html
    * marker_size= Comma-delimited list numeric marker sizes. See https://root.cern.ch/doc/master/classTAttMarker.html

    * margin_left= Double value representing the fraction of the plot size composed of a left margin
    * margin_right= Double value representing the fraction of the plot size composed of a right margin
    * margin_top= Double value representing the fraction of the plot size composed of a top margin
    * margin_bottom= Double value representing the fraction of the plot size composed of a bottom margin

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
    * Plot settings are read in from input/plot_surface.cfg. Please refer to input/template_plot_surface.cfg for full list of customizable settings.
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
        * Default: nns_response.csv

    * ***ICRP conversion factors*** - Factors to convert neutron flux [n cm^-2 s^-1] to neutron ambient dose equivalent rate [pSv s^-1]. Values were linearly interpolated from tabulated values in ICRP 74 Table A.42 (1st data column).
        * Units: pSv cm^2
        * Default: icrp_conversions.crp