#include "custom_classes.h"

#include <stdlib.h>
#include <string>

UnfoldingSettings::UnfoldingSettings() {
    norm = 1.21;
    error = 0.0;
    f_factor = 7.0;
    cutoff = 10000;
    num_poisson_samples = 50;
    // MAP specific
    beta = 0.0;
    prior = "mrp";
    // Optimize specific
    min_num_iterations = 500;
    max_num_iterations = 10000;
    iteration_increment = 50;
    min_beta = 0.0;
    max_beta = 1.0;
    parameter_of_interest = "fluence";
}

// Setter functions
void UnfoldingSettings::set_norm(double norm) {
    this->norm = norm;
}
void UnfoldingSettings::set_error(double error) {
    this->error = error;
}
void UnfoldingSettings::set_f_factor(double f_factor) {
    this->f_factor = f_factor;
}
void UnfoldingSettings::set_cutoff(int cutoff) {
    this->cutoff = cutoff;
}
void UnfoldingSettings::set_num_poisson_samples(int num_poisson_samples) {
    this->num_poisson_samples = num_poisson_samples;
}
void UnfoldingSettings::set_beta(double beta) {
    this->beta = beta;
}
void UnfoldingSettings::set_prior(std::string prior) {
    this->prior = prior;
}
void UnfoldingSettings::set_min_num_iterations(int min_num_iterations) {
    this->min_num_iterations = min_num_iterations;
}
void UnfoldingSettings::set_max_num_iterations(int max_num_iterations) {
    this->max_num_iterations = max_num_iterations;
}
void UnfoldingSettings::set_iteration_increment(int iteration_increment) {
    this->iteration_increment = iteration_increment;
}
void UnfoldingSettings::set_min_beta(double min_beta) {
    this->min_beta = min_beta;
}
void UnfoldingSettings::set_max_beta(double max_beta) {
    this->max_beta = max_beta;
}
void UnfoldingSettings::set_parameter_of_interest(std::string parameter_of_interest) {
    this->parameter_of_interest = parameter_of_interest;
}