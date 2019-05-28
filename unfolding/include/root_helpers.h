#ifndef ROOTHELPERS_H
#define ROOTHELPERS_H

#include <stdlib.h>
#include <fstream>
#include <vector>
#include <algorithm>

int plotSpectrum(std::string path_figure, 
    std::string irradiation_conditions, int num_measurements, int num_bins, 
    std::vector<double> &energy_bins, std::vector<double> &spectrum, 
    std::vector<double> &spectrum_uncertainty_upper, std::vector<double> &spectrum_uncertainty_lower
);

int plotTwoSpectra(std::string figure_file_pre, std::string figure_file_suf, 
    std::string irradiation_conditions, int num_measurements, int num_bins, 
    std::vector<double> &energy_bins, std::vector<double> &spectrum, 
    std::vector<double> &spectrum_uncertainty, std::vector<double> &second_spectrum
);

#endif