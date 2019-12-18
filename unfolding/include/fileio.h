#ifndef FILEIO_H
#define FILEIO_H

#include <stdlib.h>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include "custom_classes.h"

bool is_empty(std::ifstream& pFile);

int setSettings(std::string config_file, UnfoldingSettings &settings);

int setSpectraSettings(std::string config_file, SpectraSettings &settings);

int setPlotSettings(std::string config_file, PlotSettings &settings);

int setSurfaceSettings(std::string config_file, SurfaceSettings &settings);

bool checkStringVector(std::string item, std::vector<std::string>& allowed_items);

bool checkStringMap(std::string test_key, std::map<std::string, std::string>& test_map);

std::vector<double> getMeasurements(UnfoldingSettings &settings);

int saveSpectrumAsRow(std::string spectrum_file, int num_bins, std::string irradiation_conditions, 
    std::vector<double>& spectrum, std::vector<double> &error_lower, std::vector<double> &error_upper,
    std::vector<double>& energy_bins
);

int saveSpectrumAsColumn(std::string spectrum_file, std::string irradiation_conditions, 
    std::vector<double>& spectrum, std::vector<double>& spectrum_uncertainty, std::vector<double>& energy_bins
);

int readInputFile1D(std::string file_name, std::vector<double>& input_vector);

int readInputFile2D(std::string file_name, std::vector<std::vector<double>>& input_vector);

int readSpectra(std::string file_name, std::vector<std::string>& header_vector, std::vector<double>& energy_bins, 
    std::vector<std::vector<double>>& spectra_vector, std::vector<std::vector<double>>& error_lower_vector, 
    std::vector<std::vector<double>>& error_upper_vector, bool plot_per_mu, std::vector<int>& number_mu, 
    std::vector<int>& duration, int rows_per_spectrum
); 

int checkDimensions(int reference_size, std::string reference_string, int test_size, std::string test_string);

void stringToSVector(std::string test_string, std::vector<std::string>& result_vector);

void stringToIVector(std::string test_string, std::vector<int>& result_vector);

void stringToDVector(std::string test_string, std::vector<float>& result_vector);

int readXYYCSV(std::string file_name, std::vector<std::string>& header_vector,
    std::vector<std::vector<double>>& x_data, std::vector<std::vector<double>>& y_data
);

int readXYXYCSV(std::string file_name, std::vector<std::string>& header_vector, 
    std::vector<std::vector<double>>& x_data, std::vector<std::vector<double>>& y_data
); 

#endif