# **Neutron Spectrometry** #

This program is used to unfold data measured using a Nested Neutron Spectrometer (NNS), produce the corresponding neutron flux spectrum, and estimate the ambient neutron dose equivalent.

### Features ###

* Written in C++
* Utilizes the MLEM or MAP algorithm for unfolding
* Uncertainty estimation using poisson sampling
* Generates a report file that describes all relevant inputs/outputs of the program each time it is run
* Appends the unfolded spectrum and calculated ambient dose equivalent to output CSV files (thereby compiling subsequent results)
* Generates a plot of the unfolded spectrum and its uncertainty
* Additional functions available for data visualization, for example plotting multiple spectra in a single figure

### Setup Instructions ###

1. Clone this repository
2. Install dependencies (see dependencies.txt)
3. Compile program using make
4. Adjust configuration file (settings) if necessary
5. Provide NNS measurements
6. Run a program (including optional arguments)

For more detailed guidelines on running the program, please refer to the supplied instructions file (unfolding/README.md)