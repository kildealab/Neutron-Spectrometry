# Neutron Spectrometry

![Logo](https://github.com/McGillMedPhys/Neutron-Spectrometry/blob/master/repository_logo_figure.png)

This repository contains command-line applications to unfold data that was measured using a Nested Neutron Spectrometer (NNS), to generate plots of the resulting spectra, and other associated tasks.
<!-- Further details and instructions for running each application available [here](unfolding/instructions/). -->

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3610493.svg)](https://doi.org/10.5281/zenodo.3610493)

## Table of Contents

* [Authors](#authors)
* [Features](#features)
* [Description](#description)
* [Dependencies](#dependencies)
* [Installation](#installation)
* [List of applications](#list-of-applications)
* [Instructions](#instructions)
* [License](#license)
* [Supporting literature](#supporting-literature)

## Authors

Logan Montgomery, Georges Al Makdessi, Robert Maglieri, Felix Mathew, John Kildea

Contact email: logan.montgomery@mail.mcgill.ca

## Features

* C++
* Open source
* Command-line applications
* generate elegant plots using [ROOT](https://root.cern.ch)
* Highly customizable (using settings files)

## Description

* The primary application (`unfold_spectrum.exe`) is used to unfold NNS data via the Maximum-Likelihood Expectation&ndash;Maximization (MLEM) algorithm.
    * Multiple stopping criteria are available.
    * Our published **modified MLEM-STOP** method is included; for details see our paper: [https://doi.org/10.1016/j.nima.2020.163400](https://doi.org/10.1016/j.nima.2020.163400).
    * Spectral uncertainty estimation using a Poisson or Gaussian sampling approach.
    * Output the spectral data, a plot of the spectrum, the ambient dose equivalent, and a summary text report.
    * Easily extended to unfolding of neutron spectral data from other spectrometers (e.g. a Bonner Sphere Spectrometer).
* Additional applications to plot multiple spectra, track unfolding parameter changes with iteration number, and plot 2D data.

## Dependencies

1. C++ compiler (e.g. [GCC](https://gcc.gnu.org/))
    * Must be C++11 compatible (e.g. GCC 4.8.1)
2. GNU make ([link](https://www.gnu.org/software/make/))
3. ROOT Data Analysis Framework ([link](https://root.cern.ch/))
    * Can be installed on OSX using [Homebrew](https://brew.sh/)

**Note**: These applications were developed on OSX Mojave 10.14.6 and have been tested in Ubuntu 18.04 LTS.

## Installation

1. Download the latest version from the [releases page](https://github.com/McGillMedPhys/Neutron-Spectrometry/releases).
2. Install the [dependencies](#dependencies).
3. Execute miscellaneous initialization functions:
```
./initialization.sh
``` 
4. Compile the applications:
```
cd unfolding
make
```

## List of applications

| Application | Description |
| ----------- | ----------- |
| [`unfold_spectrum.exe`](unfolding/instructions/instructions_unfold_spectrum.md) | Read-in measured spectrometer data and unfold the neutron fluence spectrum. |
| [`plot_spectra.exe`](unfolding/instructions/instructions_plot_spectra.md) | Generate plot of one or more neutron fluence spectra. |
| [`unfold_trend.exe`](unfolding/instructions/instructions_unfold_trend.md) | Output values for a parameter of interest at each MLEM iteration. |
| [`plot_lines.exe`](unfolding/instructions/instructions_plot_lines.md) | Generate plot of one or more arbitrary sets of XY data. |

## Instructions

The general process for running each application is the same (example instructions are provided for *unfold_spectrum.exe*):

<!-- 1. Compile the applications: `make` -->
1. Enter desired settings for the application by editing the corresponding settings file (e.g. edit `input/unfold_spectrum.cfg`)
2. Prepare the data input file (e.g. edit `input/measurements.txt`)
3. Run the application (e.g. `./unfold_spectrum.exe`)

**Please click the links in the table above to view application-specific instructions** and explanations of their associated settings.

## License

This project is provided under the MIT license. See the [LICENSE file](LICENSE) for more info.

## Supporting literature

* L. Montgomery, A. Landry, F. Mathew, G. Al Makdessi, J. Kildea, "A novel MLEM stopping criterion for unfolding neutron fluence spectra in radiation therapy," *Nucl Instrum Meth A*. 2020. [https://doi.org/10.1016/j.nima.2020.163400](https://doi.org/10.1016/j.nima.2020.163400).
* L. Montgomery, M. Evans, L. Liang, R. Maglieri, J. Kildea, "The effect of the flattening filter on photoneutron production at 10 MV in the Varian TrueBeam linear accelerator," *Med Phys*. 2018. [http://dx.doi.org/10.1002/mp.13148](http://dx.doi.org/10.1002/mp.13148).
* R. Maglieri, A. Licea, M. Evans, J. Seuntjens, J. Kildea, "Measuring neutron spectra in radiotherapy using the nested neutrono spectrometer," *Med Phys*. 2015. [https://dx.doi.org/10.1118/1.4931963](https://dx.doi.org/10.1118/1.4931963).
