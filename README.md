# Neutron Spectrometry

![Logo](https://github.com/McGillMedPhys/Neutron-Spectrometry/blob/master/repository_logo_figure.png)

This repository contains command-line applications to unfold data that was measured using a Nested Neutron Spectrometer (NNS), to generate plots of the resulting spectra, and other associated tasks. Further details and instructions for running each application available [here](unfolding/README.md).

## Table of Contents

* [Features](#features)
* [Description](#description)
* [Installation](#installation)
* [Licence](#licence)
* [Supporting literature](#supporting-literature)

## Features

* C++
* Open source
* Command-line applications
* Generates elegant plots using [ROOT](https://root.cern.ch/license)
* Highly customizable (using settings files)

## Description

* Primary application is used to unfold NNS data via the Maximum-Likelihood Expectation&ndash;Maximization (MLEM) algorithm.
* Multiple stopping criteria are available.
* Our published modified MLEM-STOP method is included; for details see our peer-reviewed paper: [https://doi.org/10.1016/j.nima.2020.163400](https://doi.org/10.1016/j.nima.2020.163400).
* Spectral uncertainty estimation using a Possoin or Gaussian sampling approach.
* Output the spectral data, ambient dose equivalent, a text report summarizing the unfolding, and a plot of the spectrum.
* Complementary applications to generate elegant spectral plots, track parameter changes with iteration number, and plot 2D data.

## Installation

1. Download the latest version from the [releases page](https://github.com/McGillMedPhys/Neutron-Spectrometry/releases).
2. Install all the dependencies listed [here](unfolding/dependencies.txt).
3. Compile all the applications using make:
```
cd unfolding
make
```
4. Instructions for running specific applications are provided [here](unfolding/README.md). 

## Licence

Link will go here.

## Supporting literature

* L. Montgomery, A. Landry, F. Mathew, G. Al Makdessi, J. Kildea, "A novel MLEM stopping criterion for unfolding neutron fluence spectra in radiation therapy," *Nucl Instrum Meth A*. 2020. [https://doi.org/10.1016/j.nima.2020.163400](https://doi.org/10.1016/j.nima.2020.163400).
* L. Montgomery, M. Evans, L. Liang, R. Maglieri, J. Kildea, "The effect of the flattening filter on photoneutron production at 10 MV in the Varian TrueBeam linear accelerator," *Med Phys*. 2018. [http://dx.doi.org/10.1002/mp.13148](http://dx.doi.org/10.1002/mp.13148).
* R. Maglieri, A. Licea, M. Evans, J. Seuntjens, J. Kildea, "Measuring neutron spectra in radiotherapy using the nested neutrono spectrometer," *Med Phys*. 2015. [https://dx.doi.org/10.1118/1.4931963](https://dx.doi.org/10.1118/1.4931963).
