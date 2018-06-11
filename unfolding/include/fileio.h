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
std::vector<double> getMeasurements(std::string input_file, std::string &irradiation_conditions, double &dose_mu, double &doserate_mu, int &duration);
int saveDose(std::string dose_file, std::string irradiation_conditions, double dose, double dose_uncertainty);
int saveSpectrumAsRow(std::string spectrum_file, int num_bins, std::string irradiation_conditions, std::vector<double>& spectrum, std::vector<double>& spectrum_uncertainty, std::vector<double>& energy_bins);
int saveSpectrumAsColumn(std::string spectrum_file, std::string irradiation_conditions, std::vector<double>& spectrum, std::vector<double>& spectrum_uncertainty, std::vector<double>& energy_bins);
int readInputFile1D(std::string file_name, std::vector<double>& input_vector);
int readInputFile2D(std::string file_name, std::vector<std::vector<double>>& input_vector);
int readSpectra(std::string file_name, std::vector<std::string>& header_vector, std::vector<double>& energy_bins, std::vector<std::vector<double>>& spectra_vector, std::vector<std::vector<double>>& error_vector, bool plot_per_mu, std::vector<int>& number_mu, std::vector<int>& duration);
int checkDimensions(int reference_size, std::string reference_string, int test_size, std::string test_string);
int prepareReport(std::string report_file, std::string irradiation_conditions, std::vector<std::string> &input_files, std::vector<std::string> &input_file_flags, std::string algorithm_name, int cutoff, double error, double norm, double f_factor, double beta, int num_measurements, int num_bins, int num_poisson_samples, std::vector<double>& measurements_nc, double dose_mu, double doserate_mu, int duration, std::vector<double>& energy_bins, std::vector<double>& initial_spectrum, std::vector<std::vector<double>>& nns_response, int num_iterations, std::vector<double>& mlem_ratio, double dose, double s_dose, double total_charge, double total_flux, double total_flux_uncertainty, double avg_energy, double avg_energy_uncertainty, std::vector<double>& spectrum, std::vector<double>& spectrum_uncertainty, std::vector<double>& icrp_factors, std::string git_commit);

void stringToSVector(std::string test_string, std::vector<std::string>& result_vector);
void stringToIVector(std::string test_string, std::vector<int>& result_vector);
void stringToDVector(std::string test_string, std::vector<float>& result_vector);
int readXYYCSV(std::string file_name, std::vector<std::string>& header_vector, std::vector<double>& x_data, std::vector<std::vector<double>>& y_data);


#endif