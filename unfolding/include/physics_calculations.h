#ifndef PHYSICS_CALCULATIONS_H
#define PHYSICS_CALCULATIONS_H

#include <stdlib.h>
#include <fstream>
#include <vector>
#include <algorithm>

int processMeasurements(int num_measurements, int num_meas_per_shell, std::vector<double>& measurements, 
    std::vector<double>& std_errors);

double getMeanValueD(std::vector<double>& data);

double getSampleMeanStandardErrorD(std::vector<double>& data, double mean);

std::vector<double> normalizeResponse(int num_bins, int num_measurements, std::vector<std::vector<double>>& system_response);

std::vector<double> normalizeVector(std::vector<double>& unnormalized_vector);

double poisson(double lambda);

int runMLEM(int cutoff, double error, int num_measurements, int num_bins, std::vector<double> &measurements, 
    std::vector<double> &spectrum, std::vector<std::vector<double>>& nns_response, 
    std::vector<double> &normalized_response, std::vector<double> &mlem_ratio, 
    std::vector<double> &mlem_correction, std::vector<double> &mlem_estimate
);

int runMLEMSTOP(int cutoff, int num_measurements, int num_bins, std::vector<double> &measurements, 
    std::vector<double> &spectrum, std::vector<std::vector<double>>& nns_response, 
    std::vector<double> &normalized_response, std::vector<double> &mlem_ratio, 
    std::vector<double> &mlem_correction, std::vector<double> &mlem_estimate, double j_threshold,
    double& j_factor
);

double determineJThreshold(int num_measurements, std::vector<double>& measurements, double cps_crossover);

int runMAP(std::vector<double> &energy_correction, double beta, std::string prior, int cutoff, 
    double error, int num_measurements, int num_bins, std::vector<double> &measurements, 
    std::vector<double> &spectrum, std::vector<std::vector<double>>& nns_response, 
    std::vector<double> &normalized_response, std::vector<double> &mlem_ratio
);

double calculateDose(int num_bins, std::vector<double> &spectrum, std::vector<double> &icrp_factors);

double calculateTotalCharge(int num_measurements, std::vector<double> measurements_nc);

double calculateTotalFlux(int num_bins, std::vector<double> &spectrum);

double calculateTotalEnergyCorrection(std::vector<double> &energy_correction);

double calculateMaxRatio(int num_measurements, std::vector<double> &mlem_ratio);

double calculateAvgRatio(int num_measurements, std::vector<double> &mlem_ratio);

double calculateAverageEnergy(int num_bins, std::vector<double> &spectrum, std::vector<double> &energy_bins);

double calculateSourceStrength(int num_bins, std::vector<double> &spectrum, int duration, double dose_mu);

double calculateRMSEstimator(int size, std::vector<double> &true_vector, std::vector<double> &estimate_vector);

double calculateNRMSD(int size, std::vector<double> &true_vector, std::vector<double> &estimate_vector);

double calculateChiSquaredG(int size, std::vector<double> &true_vector, std::vector<double> &estimate_vector);

int calculateRMSD_vector(int num_samples, std::vector<double> &true_vector, 
    std::vector<std::vector<double>> &sampled_vectors, std::vector<double> &rms_differences
);

double calculateRMSD(int num_samples, double true_value, std::vector<double> &sample_vector);

double calculateSumUncertainty(int num_values, std::vector<double> &value_uncertainties);

double calculateEnergyUncertainty(int num_bins, std::vector<double> energy_bins, std::vector<double> spectrum, 
    std::vector<double> spectrum_uncertainty, double total_flux, double total_flux_uncertainty
);

double calculateJFactor(int num_measurements, std::vector<double> &measurements,
    std::vector<double> &mlem_estimate
);

double calculateJFactor2(int num_measurements, std::vector<double> &measurements,
    std::vector<double> &mlem_estimate
); 

double calculateNoise(int start_bin, int end_bin, std::vector<double>& spectrum);

double calculateChiSquared(int i_num, int num_bins, int num_measurements, std::vector<double> &spectrum, 
    std::vector<double> &measurements, std::vector<double> &mlem_ratio
);

void calculateDerivatives(std::vector<double> &derivatives, int num_points, std::vector<int> &x_data, 
    std::vector<double> &y_data
);

std::vector<double> linearSpacedDoubleVector(double a, double b, std::size_t N);

std::vector<int> linearSpacedIntegerVector(int a, int b, std::size_t N);

#endif