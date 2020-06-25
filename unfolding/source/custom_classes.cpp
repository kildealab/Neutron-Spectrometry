#include "custom_classes.h"
#include "fileio.h"
#include "physics_calculations.h"

#include <stdlib.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>

//--------------------------------------------------------------------------------------------------
// Default constructor for UnfoldingSettings
//--------------------------------------------------------------------------------------------------
UnfoldingSettings::UnfoldingSettings() {
    norm = 1.14;
    error = 0.0;
    f_factor = 7.2;
    cutoff = 15000;
    uncertainty_type = "poisson";
    num_uncertainty_samples = 50;
    num_meas_per_shell = 1;
    meas_units = "nc";
    // Measurement specs
    dose_mu = 0;
    doserate_mu = 0;
    duration = 0;
    irradiation_conditions = "radiation_details";
    // MAP specific
    beta = 0.0;
    prior = "mrp";
    // MLEM-STOP specific
    cps_crossover = 30000;
    sigma_j=0.5;
    // Optimize specific
    iteration_min = 100;
    iteration_max = 10000;
    iteration_increment = 100;
    beta_min = 1E-10;
    beta_max = 1E-8;
    parameter_of_interest = "total_fluence";
    algorithm = "mlem";
    trend_type = "cps";
    path_output_spectra = "output/output_spectra.csv";
    generate_report = 1;
    path_report = "";
    generate_figure = 1;
    path_figure = "";
    path_output_trend = "output/output_trend.csv";
    derivatives = 0;
    path_measurements = "input/measurements.txt";
    path_input_spectrum = "input/spectrum_step.csv";
    path_energy_bins = "input/energy_bins.csv";
    path_system_response = "input/response_nns_he3.csv";
    path_icrp_factors = "input/icrp_conversion_coefficients.csv";
    path_ref_spectrum = "";
}

// Apply a value to a setting:
void UnfoldingSettings::set_setting(std::string settings_name, std::string settings_value) {
    if (settings_value == ""){
        // If no setting value provided, do not apply anything
    }
    else if (settings_name == "nns_normalization")
        this->set_norm(atof(settings_value.c_str()));
    else if (settings_name == "mlem_max_error")
        this->set_error(atof(settings_value.c_str()));
    else if (settings_name == "f_factor")
        this->set_f_factor(atof(settings_value.c_str()));
    else if (settings_name == "mlem_cutoff")
        this->set_cutoff(atoi(settings_value.c_str()));
    else if (settings_name == "uncertainty_type")
        this->set_uncertainty_type(settings_value);
    else if (settings_name == "num_uncertainty_samples")
        this->set_num_uncertainty_samples(atoi(settings_value.c_str()));
    else if (settings_name == "num_meas_per_shell")
        this->set_num_meas_per_shell(atoi(settings_value.c_str()));
    else if (settings_name == "meas_units")
        this->set_meas_units(settings_value);
    else if (settings_name == "dose_mu")
        this->set_dose_mu(atoi(settings_value.c_str()));
    else if (settings_name == "doserate_mu")
        this->set_doserate_mu(atoi(settings_value.c_str()));
    else if (settings_name == "duration")
        this->set_duration(atoi(settings_value.c_str()));
    else if (settings_name == "irradiation_conditions")
        this->set_irradiation_conditions(settings_value);
    else if (settings_name == "beta")
        this->set_beta(atof(settings_value.c_str()));
    else if (settings_name == "prior")
        this->set_prior(settings_value);
    else if (settings_name == "cps_crossover")
        this->set_cps_crossover(atoi(settings_value.c_str()));
    else if (settings_name == "sigma_j")
        this->set_sigma_j(atof(settings_value.c_str()));
    else if (settings_name == "iteration_min")
        this->set_iteration_min(atoi(settings_value.c_str()));
    else if (settings_name == "iteration_max")
        this->set_iteration_max(atoi(settings_value.c_str()));
    else if (settings_name == "iteration_increment")
        this->set_iteration_increment(atoi(settings_value.c_str()));
    else if (settings_name == "beta_min")
        this->set_beta_min(atof(settings_value.c_str()));
    else if (settings_name == "beta_max")
        this->set_beta_max(atof(settings_value.c_str()));
    else if (settings_name == "parameter_of_interest")
        this->set_parameter_of_interest(settings_value);
    else if (settings_name == "algorithm")
        this->set_algorithm(settings_value);
    else if (settings_name == "trend_type")
        this->set_trend_type(settings_value);
    else if (settings_name == "path_output_spectra")
        this->set_path_output_spectra(settings_value);
    else if (settings_name == "generate_report")
        this->set_generate_report(atoi(settings_value.c_str()));
    else if (settings_name == "path_report")
        this->set_path_report(settings_value);
    else if (settings_name == "generate_figure")
        this->set_generate_figure(atoi(settings_value.c_str()));
    else if (settings_name == "path_figure")
        this->set_path_figure(settings_value);
    else if (settings_name == "path_output_trend")
        this->set_path_output_trend(settings_value);
    else if (settings_name == "derivatives")
        this->set_derivatives(atoi(settings_value.c_str()));
    else if (settings_name == "path_measurements")
        this->set_path_measurements(settings_value);
    else if (settings_name == "path_input_spectrum")
        this->set_path_input_spectrum(settings_value);
    else if (settings_name == "path_energy_bins")
        this->set_path_energy_bins(settings_value);
    else if (settings_name == "path_system_response")
        this->set_path_system_response(settings_value);
    else if (settings_name == "path_icrp_factors")
        this->set_path_icrp_factors(settings_value);
    else if (settings_name == "path_ref_spectrum")
        this->set_path_ref_spectrum(settings_value);
    else
        throw std::logic_error("Unrecognized setting: " + settings_name 
            + ". Please refer to the README for allowed settings");

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
void UnfoldingSettings::set_uncertainty_type(std::string uncertainty_type) {
    this->uncertainty_type = uncertainty_type;
}
void UnfoldingSettings::set_num_uncertainty_samples(int num_uncertainty_samples) {
    this->num_uncertainty_samples = num_uncertainty_samples;
}
void UnfoldingSettings::set_num_meas_per_shell(int num_meas_per_shell) {
    this->num_meas_per_shell = num_meas_per_shell;
}
void UnfoldingSettings::set_meas_units(std::string meas_units) {
    this->meas_units = meas_units;
}
void UnfoldingSettings::set_dose_mu(int dose_mu) {
    this->dose_mu = dose_mu;
}
void UnfoldingSettings::set_doserate_mu(int doserate_mu) {
    this->doserate_mu = doserate_mu;
}
void UnfoldingSettings::set_duration(int duration) {
    this->duration = duration;
}
void UnfoldingSettings::set_irradiation_conditions(std::string irradiation_conditions) {
    this->irradiation_conditions = irradiation_conditions;
}
void UnfoldingSettings::set_beta(double beta) {
    this->beta = beta;
}
void UnfoldingSettings::set_prior(std::string prior) {
    this->prior = prior;
}
void UnfoldingSettings::set_cps_crossover(int cps_crossover) {
    this->cps_crossover = cps_crossover;
}
void UnfoldingSettings::set_sigma_j(double sigma_j) {
    this->sigma_j = sigma_j;
}
void UnfoldingSettings::set_iteration_min(int iteration_min) {
    this->iteration_min = iteration_min;
}
void UnfoldingSettings::set_iteration_max(int iteration_max) {
    this->iteration_max = iteration_max;
}
void UnfoldingSettings::set_iteration_increment(int iteration_increment) {
    this->iteration_increment = iteration_increment;
}
void UnfoldingSettings::set_beta_min(double beta_min) {
    this->beta_min = beta_min;
}
void UnfoldingSettings::set_beta_max(double beta_max) {
    this->beta_max = beta_max;
}
void UnfoldingSettings::set_parameter_of_interest(std::string parameter_of_interest) {
    this->parameter_of_interest = parameter_of_interest;
}
void UnfoldingSettings::set_algorithm(std::string algorithm) {
    this->algorithm = algorithm;
}
void UnfoldingSettings::set_trend_type(std::string trend_type) {
    this->trend_type = trend_type;
}
void UnfoldingSettings::set_path_output_spectra(std::string path_output_spectra) {
    this->path_output_spectra = path_output_spectra;
}
void UnfoldingSettings::set_generate_report(int generate_report) {
    this->generate_report = generate_report;
}
void UnfoldingSettings::set_path_report(std::string path_report) {
    this->path_report = path_report;
}
void UnfoldingSettings::set_generate_figure(int generate_figure) {
    this->generate_figure = generate_figure;
}
void UnfoldingSettings::set_path_figure(std::string path_figure) {
    this->path_figure = path_figure;
}
void UnfoldingSettings::set_path_output_trend(std::string path_output_trend) {
    this->path_output_trend = path_output_trend;
}
void UnfoldingSettings::set_derivatives(int derivatives) {
    this->derivatives = derivatives;
}
void UnfoldingSettings::set_path_measurements(std::string path_measurements) {
    this->path_measurements = path_measurements;
}
void UnfoldingSettings::set_path_input_spectrum(std::string path_input_spectrum) {
    this->path_input_spectrum = path_input_spectrum;
}
void UnfoldingSettings::set_path_energy_bins(std::string path_energy_bins) {
    this->path_energy_bins = path_energy_bins;
}
void UnfoldingSettings::set_path_system_response(std::string path_system_response) {
    this->path_system_response = path_system_response;
}
void UnfoldingSettings::set_path_icrp_factors(std::string path_icrp_factors) {
    this->path_icrp_factors = path_icrp_factors;
}
void UnfoldingSettings::set_path_ref_spectrum(std::string path_ref_spectrum) {
    this->path_ref_spectrum = path_ref_spectrum;
}


//--------------------------------------------------------------------------------------------------
// Default constructor for UnfoldingSettings
//--------------------------------------------------------------------------------------------------
UnfoldingReport::UnfoldingReport() {
    path = "output/report.txt";
}

void UnfoldingReport::set_path(std::string path) {
    this->path = path;
}
void UnfoldingReport::set_irradiation_conditions(std::string irradiation_conditions) {
    this->irradiation_conditions = irradiation_conditions;
}
void UnfoldingReport::set_input_files(std::vector<std::string>& input_files) {
    this->input_files = input_files;
}
void UnfoldingReport::set_input_file_flags(std::vector<std::string>& input_file_flags) {
    this->input_file_flags = input_file_flags;
}
void UnfoldingReport::set_cutoff(int cutoff) {
    this->cutoff = cutoff;
}
void UnfoldingReport::set_error(double error) {
    this->error = error;
}
void UnfoldingReport::set_norm(double norm) {
    this->norm = norm;
}
void UnfoldingReport::set_f_factor(double f_factor) {
    this->f_factor = f_factor;
}
void UnfoldingReport::set_num_measurements(int num_measurements) {
    this->num_measurements = num_measurements;
}
void UnfoldingReport::set_num_bins(int num_bins) {
    this->num_bins = num_bins;
}
void UnfoldingReport::set_uncertainty_type(std::string uncertainty_type) {
    this->uncertainty_type = uncertainty_type;
}
void UnfoldingReport::set_num_uncertainty_samples(int num_uncertainty_samples) {
    this->num_uncertainty_samples = num_uncertainty_samples;
}
void UnfoldingReport::set_git_commit(std::string git_commit) {
    this->git_commit = git_commit;
}
void UnfoldingReport::set_measurements(std::vector<double>& measurements) {
    this->measurements = measurements;
}
void UnfoldingReport::set_measurements_nc(std::vector<double>& measurements_nc) {
    this->measurements_nc = measurements_nc;
}
void UnfoldingReport::set_dose_mu(double dose_mu) {
    this->dose_mu = dose_mu;
}
void UnfoldingReport::set_doserate_mu(double doserate_mu) {
    this->doserate_mu = doserate_mu;
}
void UnfoldingReport::set_duration(int duration) {
    this->duration = duration;
}
void UnfoldingReport::set_meas_units(std::string meas_units) {
    this->meas_units = meas_units;
}
void UnfoldingReport::set_initial_spectrum(std::vector<double>& initial_spectrum) {
    this->initial_spectrum = initial_spectrum;
}
void UnfoldingReport::set_energy_bins(std::vector<double>& energy_bins) {
    this->energy_bins = energy_bins;
}
void UnfoldingReport::set_nns_response(std::vector<std::vector<double>>& nns_response) {
    this->nns_response = nns_response;
}
void UnfoldingReport::set_icrp_factors(std::vector<double>& icrp_factors) {
    this->icrp_factors = icrp_factors;
}
void UnfoldingReport::set_spectrum(std::vector<double>& spectrum) {
    this->spectrum = spectrum;
}
void UnfoldingReport::set_spectrum_uncertainty_upper(std::vector<double>& spectrum_uncertainty_upper) {
    this->spectrum_uncertainty_upper = spectrum_uncertainty_upper;
}
void UnfoldingReport::set_spectrum_uncertainty_lower(std::vector<double>& spectrum_uncertainty_lower) {
    this->spectrum_uncertainty_lower = spectrum_uncertainty_lower;
}
void UnfoldingReport::set_num_iterations(int num_iterations) {
    this->num_iterations = num_iterations;
}
void UnfoldingReport::set_mlem_ratio(std::vector<double>& mlem_ratio) {
    this->mlem_ratio = mlem_ratio;
}
void UnfoldingReport::set_dose(double dose) {
    this->dose = dose;
}
void UnfoldingReport::set_dose_uncertainty_upper(double dose_uncertainty_upper) {
    this->dose_uncertainty_upper = dose_uncertainty_upper;
}
void UnfoldingReport::set_dose_uncertainty_lower(double dose_uncertainty_lower) {
    this->dose_uncertainty_lower = dose_uncertainty_lower;
}
void UnfoldingReport::set_total_flux(double total_flux) {
    this->total_flux = total_flux;
}
void UnfoldingReport::set_total_flux_uncertainty_upper(double total_flux_uncertainty_upper) {
    this->total_flux_uncertainty_upper = total_flux_uncertainty_upper;
}
void UnfoldingReport::set_total_flux_uncertainty_lower(double total_flux_uncertainty_lower) {
    this->total_flux_uncertainty_lower = total_flux_uncertainty_lower;
}
void UnfoldingReport::set_avg_energy(double avg_energy) {
    this->avg_energy = avg_energy;
}
void UnfoldingReport::set_avg_energy_uncertainty_upper(double avg_energy_uncertainty_upper) {
    this->avg_energy_uncertainty_upper = avg_energy_uncertainty_upper;
}
void UnfoldingReport::set_avg_energy_uncertainty_lower(double avg_energy_uncertainty_lower) {
    this->avg_energy_uncertainty_lower = avg_energy_uncertainty_lower;
}
void UnfoldingReport::set_algorithm(std::string algorithm) {
    this->algorithm = algorithm;
}
void UnfoldingReport::set_cps_crossover(int cps_crossover) {
    this->cps_crossover = cps_crossover;
}
void UnfoldingReport::set_j_threshold(double j_threshold) {
    this->j_threshold = j_threshold;
}
void UnfoldingReport::set_j_final(double j_final) {
    this->j_final = j_final;
}
void UnfoldingReport::set_j_manager_low(UncertaintyManagerJ j_manager_low) {
    this->j_manager_low = j_manager_low;
}
void UnfoldingReport::set_j_manager_high(UncertaintyManagerJ j_manager_high) {
    this->j_manager_high = j_manager_high;
}
void UnfoldingReport::set_num_toss(int num_toss) {
    this->num_toss = num_toss;
}

//----------------------------------------------------------------------------------------------
// Prepare summary report of unfolding
//----------------------------------------------------------------------------------------------
void UnfoldingReport::prepare_report() {
    std::ofstream rfile(path);

    report_header(rfile);
    report_settings(rfile);
    report_measurement_info(rfile);
    report_inputs(rfile);
    report_mlem_info(rfile);
    report_results(rfile);

    rfile.close();
}

//----------------------------------------------------------------------------------------------
// Header
//----------------------------------------------------------------------------------------------
void UnfoldingReport::report_header(std::ofstream& rfile) {
    rfile << HEADER_DIVIDE;
    rfile << "Neutron Unfolding Report\n\n";
    rfile << std::left << std::setw(sw) << "Irradiation Specs: " << irradiation_conditions << "\n";
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    rfile << std::left << std::setw(sw) << "Date report was generated: " 
        << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << "\n";
    if (git_commit != ""){
        rfile << std::left << std::setw(sw) << "Git commit number: " << git_commit << "\n";
    }
    rfile << "Input arguments (files) used:\n";
    for (int i=0; i<input_files.size(); i++) {
        std::string tempstring = "    " + input_file_flags[i];
        rfile << std::left <<std::setw(sw) << tempstring << input_files[i] << "\n";
    }
    rfile << HEADER_DIVIDE << "\n";
}

//----------------------------------------------------------------------------------------------
// Settings
//----------------------------------------------------------------------------------------------
void UnfoldingReport::report_settings(std::ofstream& rfile) {
    rfile << "Settings\n\n";
    rfile << std::left << std::setw(sw) << "MLEM max # of iterations:" << cutoff << "\n";
    rfile << std::left << std::setw(sw) << "MLEM target ratio:" << error << "\n";
    rfile << std::left << std::setw(sw) << "NNS normalization factor:" << norm << "\n";
    rfile << std::left << std::setw(sw) << "NNS calibration factor:" << f_factor << " fA/cps\n";
    rfile << std::left << std::setw(sw) << "Uncertainty type:" << uncertainty_type << " fA/cps\n";
    rfile << std::left << std::setw(sw) << "# of uncertainty samples:" << num_uncertainty_samples << "\n";
    if (algorithm == "mlemstop") {
        rfile << std::left << std::setw(sw) << "Crossover CPS value:" << cps_crossover << "\n";
        rfile << std::left << std::setw(sw) << "J threshold:" << j_threshold << "\n";
    }
    rfile << SECTION_DIVIDE;
}

//----------------------------------------------------------------------------------------------
// Measurement info
//----------------------------------------------------------------------------------------------
void UnfoldingReport::report_measurement_info(std::ofstream& rfile) {
    std::string units_string;
    if (meas_units == "nc")
        units_string = "Charge (nC)";
    else
        units_string = "CPS";
    
    rfile << "Measurement\n\n";
    rfile << std::left << std::setw(sw) << "Delivered dose:" << dose_mu << " MU\n";
    rfile << std::left << std::setw(sw) << "Delivered doserate:" << doserate_mu << " MU/min\n";
    rfile << std::left << std::setw(sw) << "Measurement duration:" << duration << " s\n\n";
    // rfile << "Measured Data (measurement duration: " << duration << "s)\n\n";

    if (meas_units == "nc") { 
        rfile << std::left << std::setw(cw) << "# of moderators" << std::setw(cw) << "CPS" << "Charge (nC)" "\n";
        rfile << std::left << std::setw(cw) << COLSTRING << std::setw(cw) << COLSTRING << COLSTRING << "\n";
        for (int i=0; i<num_measurements; i++) {
            rfile << std::left << std::setw(cw) << i << std::setw(cw) << round(measurements[i]) << measurements_nc[i] << "\n";
        }
    }
    else { 
        rfile << std::left << std::setw(cw) << "# of moderators" << "CPS" "\n";
        rfile << std::left << std::setw(cw) << COLSTRING << COLSTRING << "\n";
        for (int i=0; i<num_measurements; i++) {
            rfile << std::left << std::setw(cw) << i << round(measurements[i]) << "\n";
        }
    }

    rfile << SECTION_DIVIDE;
}

//----------------------------------------------------------------------------------------------
// Inputs
//----------------------------------------------------------------------------------------------
void UnfoldingReport::report_inputs(std::ofstream& rfile) {
    rfile << "Inputs (Number of energy bins: " << num_bins << ")\n\n";
    rfile << std::left << std::setw(cw) << "Energy bins" << std::setw(cw) << "Input spectrum" 
        << "| NNS Response by # of moderators (cm^2)\n";
    rfile << std::left << std::setw(cw) << "(MeV)" << std::setw(cw) << "(n cm^-2 s^-1)" << "| ";
    for (int j=0; j<num_measurements; j++) {
        rfile << std::left << std::setw(rw) << j;
    }
    rfile << "\n";
    rfile << std::left << std::setw(cw) << COLSTRING << std::setw(cw) << COLSTRING << "--";
    for (int j=0; j<num_measurements; j++) {
        rfile << "---------";
    }
    rfile << "\n";

    for (int i=0; i<num_bins; i++) {
        rfile << std::left << std::setw(cw) << energy_bins[i] << std::setw(cw) << initial_spectrum[i] << "| ";
        for (int j=0; j<num_measurements; j++) {
            rfile << std::left << std::setw(rw) << nns_response[j][i];
        }
        rfile << "\n";
    }
    rfile << SECTION_DIVIDE;
}

//----------------------------------------------------------------------------------------------
// MLEM Processing
//----------------------------------------------------------------------------------------------
void UnfoldingReport::report_mlem_info(std::ofstream& rfile) {
    rfile << "Unfolding information\n\n";
    rfile << std::left << std::setw(sw) << "Algorithm: " << algorithm << "\n";
    if (algorithm == "map") {
        rfile << std::left << std::setw(sw) << "MAP beta value: " << beta << "\n";
    }
    rfile << std::left << std::setw(sw) << "# of iterations: " << num_iterations << "/" << cutoff << "\n\n";
    if (algorithm == "mlemstop") {
        rfile << std::left << std::setw(sw) << "final J value: " << j_final << "/" << j_threshold << "\n\n";
    }
    if (algorithm == "mlemstop") {
        rfile << std::left << std::setw(sw) << "# samples tossed: " << num_toss << "/" << num_uncertainty_samples+num_toss << "\n\n";
    }
    rfile << "Final unfolding ratio = measured charge / estimated charge:\n";
    int thw = 13; // NNS response column width
    //row 1
    rfile << std::left << std::setw(thw) << "# moderators" << "| ";
    for (int j=0; j<num_measurements; j++) {
        rfile << std::left << std::setw(rw) << j;
    }
    rfile << "\n";
    // row 2
    rfile << std::left << std::setw(thw) << "-------------|-";
    for (int j=0; j<num_measurements; j++) {
        rfile << "---------";
    }
    rfile << "\n";
    // row thw
    rfile << std::left << std::setw(thw) << "ratio" << "| ";
    for (int j=0; j<num_measurements; j++) {
        rfile << std::left << std::setw(rw) << mlem_ratio[j];
    }
    rfile << "\n";
    rfile << SECTION_DIVIDE;
}

//----------------------------------------------------------------------------------------------
// Results
//----------------------------------------------------------------------------------------------
void UnfoldingReport::report_results(std::ofstream& rfile) {
    rfile << "Results\n\n";
    rfile << std::left << std::setw(sw) << "Ambient dose equivalent:" << dose << " mSv/hr\n";
    rfile << std::left << std::setw(sw) << "Upper uncertainty:" << dose_uncertainty_upper << " mSv/hr\n";
    rfile << std::left << std::setw(sw) << "Lower uncertainty:" << dose_uncertainty_lower << " mSv/hr\n\n";
    // rfile << std::left << std::setw(sw) << "Total measured charge:" << total_charge << " nC\n\n";
    if (uncertainty_type == "j_bounds") {
        rfile << std::left << std::setw(sw) << "Upper uncertainty:" << j_manager_high.j_factor << "/" 
            << j_manager_high.j_threshold << " (k=" << j_manager_high.num_iterations << ")\n";
        rfile << std::left << std::setw(sw) << "Lower uncertainty:" << j_manager_low.j_factor << "/" 
            << j_manager_low.j_threshold << " (k=" << j_manager_low.num_iterations << ")\n\n";
    }
    rfile << std::left << std::setw(sw) << "Integrated neutron flux:" << total_flux << " n cm^-2 s^-1\n";
    rfile << std::left << std::setw(sw) << "Upper uncertainty:" << total_flux_uncertainty_upper << " n cm^-2 s^-1\n";
    rfile << std::left << std::setw(sw) << "Lower uncertainty:" << total_flux_uncertainty_lower << " n cm^-2 s^-1\n\n";
    rfile << std::left << std::setw(sw) << "Average neutron energy:" << avg_energy << " MeV\n";
    rfile << std::left << std::setw(sw) << "Upper uncertainty:" << avg_energy_uncertainty_upper << " MeV\n";
    rfile << std::left << std::setw(sw) << "Lower uncertainty:" << avg_energy_uncertainty_lower << " MeV\n\n";

    rfile << std::left << std::setw(cw) << "Energy bins" << std::setw(cw) << "Unfolded spectrum" << std::setw(cw) 
        << "+ Uncertainty" << std::setw(cw) << "- Uncertainty" << std::setw(cw) << "| ICRP H factor" << "Ambient Dose Equiv.\n";
    rfile << std::left << std::setw(cw) << "(MeV)" << std::setw(cw) << "(n cm^-2 s^-1)" << std::setw(cw) 
        << "(n cm^-2 s^-1)" << std::setw(cw) << "(n cm^-2 s^-1)" << std::setw(cw) << "| (pSv/cm^2)" << "(mSv/hr)\n";;
    rfile << std::left << std::setw(cw) << COLSTRING << std::setw(cw) << COLSTRING << std::setw(cw) 
        << COLSTRING << std::setw(cw) << COLSTRING << std::setw(cw) << COLSTRING << COLSTRING << "\n";
    for (int i=0; i<num_bins; i++) {
        std::ostringstream icrp_string;
        icrp_string << "| " <<icrp_factors[i];
        double subdose = spectrum[i]*icrp_factors[i]*3600*(1e-9);
        rfile << std::left << std::setw(cw) << energy_bins[i] << std::setw(cw) << spectrum[i] << std::setw(cw) 
            << spectrum_uncertainty_upper[i] << std::setw(cw) << spectrum_uncertainty_lower[i] << std::setw(26) 
            << icrp_string.str() << subdose << "\n";
    }
}

//--------------------------------------------------------------------------------------------------
// Default constructor to be used to create J Uncertainty Manager object. Required when used as
// a class member in other classes
//--------------------------------------------------------------------------------------------------
UncertaintyManagerJ::UncertaintyManagerJ() {
    j_threshold = 0;
    j_factor = 0;
    num_iterations = 0;
    dose_uncertainty = 0;
}

//--------------------------------------------------------------------------------------------------
// Constructor to be used to create J Uncertainty Manager object. This should be used in practice
//--------------------------------------------------------------------------------------------------
UncertaintyManagerJ::UncertaintyManagerJ(double original_j_threshold, double sigma_j) {
    j_threshold = original_j_threshold*sigma_j;
    j_factor = 0;
    num_iterations = 0;
    dose_uncertainty = 0;
}

//--------------------------------------------------------------------------------------------------
// Method used to calculate a single upper or lower uncertainty on a neutron fluence spectrum. The
// MLEM-STOP method is used to determine the iteration number at which the J-threshold (scaled by
// sigma_J) is attained. The spectral difference between this spectrum and the actual MLEM-STOP
// spectrum is taken to be the uncertainty.
//--------------------------------------------------------------------------------------------------
void UncertaintyManagerJ::determineSpectrumUncertainty(std::vector<double> &mlemstop_spectrum, 
    int cutoff, int num_measurements, int num_bins, std::vector<double> &measurements, 
    std::vector<std::vector<double>> &nns_response, std::vector<double> &normalized_response,
    std::vector<double> &initial_spectrum) 
{
    this->bound_spectrum = initial_spectrum;

    std::vector<double> mlem_ratio;
    std::vector<double> mlem_correction;
    std::vector<double> mlem_estimate;

    num_iterations = runMLEMSTOP(cutoff, num_measurements, num_bins, measurements,
        this->bound_spectrum, nns_response, normalized_response, mlem_ratio, mlem_correction, mlem_estimate,
        this->j_threshold, this->j_factor
    );

    for (int i_bin = 0; i_bin < num_bins; i_bin++) {
        spectrum_uncertainty.push_back(abs(bound_spectrum[i_bin]-mlemstop_spectrum[i_bin]));
    }
}

//--------------------------------------------------------------------------------------------------
// Calculate the uncertainty in dose by calculating the dose at the bound spectrum, and subtracting
// the dose for the actual MLEM-STOP spectrum
//--------------------------------------------------------------------------------------------------
void UncertaintyManagerJ::determineDoseUncertainty(double dose, std::vector<double> &mlemstop_spectrum, int num_bins, 
    std::vector<double> &icrp_factors)
{
    double dose_bound = calculateDose(num_bins,this->bound_spectrum,icrp_factors);
    this->dose_uncertainty = abs(dose_bound-dose);
}

//--------------------------------------------------------------------------------------------------
// Default Constructor for SpectraSettings
//--------------------------------------------------------------------------------------------------
SpectraSettings::SpectraSettings() {
    path_input_data = "output/output_spectra.csv";
    path_output_figure = "output/output_spectra.png";
    title = "";
    x_label = "Energy (MeV)";
    y_label = "Fluence (n #upoint cm^{-2} s^{-1})";
    x_min = 0; // if plotting detects max and min are the same, use default limits
    x_max = 0;
    y_min = 0; // if plotting detects max and min are the same, use default limits
    y_max = 0;
    x_res = 3200;
    y_res = 2400;
    y_num_divs = 0;
    y_digits_max = 3;
    legend_entries = {};
    color_series = {"#000000","#C63822","#607FD5","#55A961","#75298e","#ce884a","#33889b"};
    color_error = {"#333333","#E79A9F","#6B8EF0","#69CF77","#ad4ace","#e2a56f","#69c6db"};
    grayscale = 0;
    show_error = {1};
    error_style = "E2";
    error_fill_style = 3001;
    rows_per_spectrum = 3;
    line_style = {1};
    line_width = {5};
    border_width = 5;
    legend = 1;
    legend_coords = {0.15,0.75,0.4,0.85};
    textbox = 0;
    textbox_coords = {0.15,0.4,0.4,0.6};
    textbox_text = {};
    plot_per_mu = 0;
    number_mu = {0};
    duration = {0};
    normalize = 0;
    margin_left = 0.1;
    margin_right = 0.1;
    margin_top = 0.1;
    margin_bottom = 0.1;
    font_size = 1.0;
    font_size_legend = 1.0;
    font_size_axis_labels = 1.0;
    font_size_axis_tick_labels = 1.0;
    font_size_title = 1.0;
    x_label_offset = 1.4;
    y_label_offset = 1.4;
}

// Apply a value to a setting:
void SpectraSettings::set_setting(std::string settings_name, std::string settings_value) {
    if (settings_value == ""){
        // If no setting value provided, do not apply anything
    }
    else if (settings_name == "path_input_data")
        this->set_path_input_data(settings_value);
    else if (settings_name == "path_output_figure")
        this->set_path_output_figure(settings_value);
    else if (settings_name == "title")
        this->set_title(settings_value);
    else if (settings_name == "x_label")
        this->set_x_label(settings_value);
    else if (settings_name == "y_label")
        this->set_y_label(settings_value);
    else if (settings_name == "x_min")
        this->set_x_min(settings_value);
    else if (settings_name == "x_max")
        this->set_x_max(settings_value);
    else if (settings_name == "y_min")
        this->set_y_min(settings_value);
    else if (settings_name == "y_max")
        this->set_y_max(settings_value);
    else if (settings_name == "x_res")
        this->set_x_res(settings_value);
    else if (settings_name == "y_res")
        this->set_y_res(settings_value);
    else if (settings_name == "y_num_divs")
        this->set_y_num_divs(settings_value);
    else if (settings_name == "y_digits_max")
        this->set_y_digits_max(settings_value);
    else if (settings_name == "legend_entries")
        this->set_legend_entries(settings_value);
    else if (settings_name == "color_series")
        this->set_color_series(settings_value);
    else if (settings_name == "color_error")
        this->set_color_error(settings_value);
    else if (settings_name == "grayscale")
        this->set_grayscale(settings_value);
    else if (settings_name == "show_error")
        this->set_show_error(settings_value);
    else if (settings_name == "error_style")
        this->set_error_style(settings_value);
    else if (settings_name == "error_fill_style")
        this->set_error_fill_style(settings_value);
    else if (settings_name == "rows_per_spectrum")
        this->set_rows_per_spectrum(settings_value);
    else if (settings_name == "line_style")
        this->set_line_style(settings_value);
    else if (settings_name == "line_width")
        this->set_line_width(settings_value);
    else if (settings_name == "border_width")
        this->set_border_width(settings_value);
    else if (settings_name == "legend")
        this->set_legend(settings_value);
    else if (settings_name == "legend_coords")
        this->set_legend_coords(settings_value);
    else if (settings_name == "textbox")
        this->set_textbox(settings_value);
    else if (settings_name == "textbox_coords")
        this->set_textbox_coords(settings_value);
    else if (settings_name == "textbox_text")
        this->set_textbox_text(settings_value);
    else if (settings_name == "plot_per_mu")
        this->set_plot_per_mu(settings_value);
    else if (settings_name == "number_mu")
        this->set_number_mu(settings_value);
    else if (settings_name == "duration")
        this->set_duration(settings_value);
    else if (settings_name == "normalize")
        this->set_normalize(settings_value);
    else if (settings_name == "margin_left")
        this->set_margin_left(settings_value);
    else if (settings_name == "margin_right")
        this->set_margin_right(settings_value);
    else if (settings_name == "margin_top")
        this->set_margin_top(settings_value);
    else if (settings_name == "margin_bottom")
        this->set_margin_bottom(settings_value);
    else if (settings_name == "font_size")
        this->set_font_size(settings_value);
    else if (settings_name == "font_size_legend")
        this->set_font_size_legend(settings_value);
    else if (settings_name == "font_size_axis_labels")
        this->set_font_size_axis_labels(settings_value);
    else if (settings_name == "font_size_axis_tick_labels")
        this->set_font_size_axis_tick_labels(settings_value);
    else if (settings_name == "font_size_title")
        this->set_font_size_title(settings_value);
    else if (settings_name == "x_label_offset")
        this->set_x_label_offset(settings_value);
    else if (settings_name == "y_label_offset")
        this->set_y_label_offset(settings_value);
    // else
    //     throw std::logic_error("Unrecognized setting: " + settings_name 
    //      + ". Please refer to the README for allowed settings");
}

void SpectraSettings::set_path_input_data(std::string path_input_data) {
    this->path_input_data = path_input_data;
}
void SpectraSettings::set_path_output_figure(std::string path_output_figure) {
    this->path_output_figure = path_output_figure;
}
void SpectraSettings::set_title(std::string title) {
    this->title = title;
}
void SpectraSettings::set_x_label(std::string x_label) {
    this->x_label = x_label;
}
void SpectraSettings::set_y_label(std::string y_label) {
    this->y_label = y_label;
}
void SpectraSettings::set_x_min(std::string x_min) {
    this->x_min = stod(x_min);
}
void SpectraSettings::set_x_max(std::string x_max) {
    this->x_max = stod(x_max);
}
void SpectraSettings::set_y_min(std::string y_min) {
    this->y_min = stod(y_min);
}
void SpectraSettings::set_y_max(std::string y_max) {
    this->y_max = stod(y_max);
}
void SpectraSettings::set_x_res(std::string x_res) {
    this->x_res = stoi(x_res);
}
void SpectraSettings::set_y_res(std::string y_res) {
    this->y_res = stoi(y_res);
}
void SpectraSettings::set_y_num_divs(std::string y_num_divs) {
    this->y_num_divs = stoi(y_num_divs);
}
void SpectraSettings::set_y_digits_max(std::string y_digits_max) {
    this->y_digits_max = stoi(y_digits_max);
}
void SpectraSettings::set_legend_entries(std::string legend_entries) {
    stringToSVector(legend_entries,this->legend_entries);
}
void SpectraSettings::set_color_series(std::string color_series) {
    stringToSVector(color_series,this->color_series);
}
void SpectraSettings::set_color_error(std::string color_error) {
    stringToSVector(color_error,this->color_error);
}
void SpectraSettings::set_grayscale(std::string grayscale) {
    this->grayscale = stoi(grayscale);
} 
void SpectraSettings::set_show_error(std::string show_error) {
    stringToIVector(show_error,this->show_error);
}
void SpectraSettings::set_error_style(std::string error_style) {
    this->error_style = error_style;
}
void SpectraSettings::set_error_fill_style(std::string error_fill_style) {
    this->error_fill_style = stoi(error_fill_style);
}
void SpectraSettings::set_rows_per_spectrum(std::string rows_per_spectrum) {
    this->rows_per_spectrum = stoi(rows_per_spectrum);
}
void SpectraSettings::set_line_style(std::string line_style) {
    stringToIVector(line_style,this->line_style);
}
void SpectraSettings::set_line_width(std::string line_width) {
    stringToIVector(line_width,this->line_width);
}
void SpectraSettings::set_border_width(std::string border_width) {
    this->border_width = stoi(border_width);
}
void SpectraSettings::set_legend(std::string legend) {
    this->legend = stoi(legend);
} 
void SpectraSettings::set_legend_coords(std::string legend_coords) {
    stringToDVector(legend_coords,this->legend_coords);
}
void SpectraSettings::set_textbox(std::string textbox) {
    this->textbox = stoi(textbox);
} 
void SpectraSettings::set_textbox_coords(std::string textbox_coords) {
    stringToDVector(textbox_coords,this->textbox_coords);
}
void SpectraSettings::set_textbox_text(std::string textbox_text) {
    stringToSVector(textbox_text,this->textbox_text);
}
void SpectraSettings::set_plot_per_mu(std::string plot_per_mu) {
    this->plot_per_mu = stoi(plot_per_mu);
} 
void SpectraSettings::set_number_mu(std::string number_mu) {
    stringToIVector(number_mu,this->number_mu);
} 
void SpectraSettings::set_duration(std::string duration) {
    stringToIVector(duration,this->duration);
} 
void SpectraSettings::set_normalize(std::string normalize) {
    this->normalize = stoi(normalize);
} 
void SpectraSettings::set_margin_left(std::string margin_left) {
    this->margin_left = stod(margin_left);
}
void SpectraSettings::set_margin_right(std::string margin_right) {
    this->margin_right = stod(margin_right);
}
void SpectraSettings::set_margin_top(std::string margin_top) {
    this->margin_top = stod(margin_top);
}
void SpectraSettings::set_margin_bottom(std::string margin_bottom) {
    this->margin_bottom = stod(margin_bottom);
}
void SpectraSettings::set_font_size(std::string font_size) {
    this->font_size = stod(font_size);
}
void SpectraSettings::set_font_size_legend(std::string font_size_legend) {
    this->font_size_legend = stod(font_size_legend);
}
void SpectraSettings::set_font_size_axis_labels(std::string font_size_axis_labels) {
    this->font_size_axis_labels = stod(font_size_axis_labels);
}
void SpectraSettings::set_font_size_axis_tick_labels(std::string font_size_axis_tick_labels) {
    this->font_size_axis_tick_labels = stod(font_size_axis_tick_labels);
}
void SpectraSettings::set_font_size_title(std::string font_size_title) {
    this->font_size_title = stod(font_size_title);
}
void SpectraSettings::set_font_size_textbox(std::string font_size_textbox) {
    this->font_size_textbox = stod(font_size_textbox);
}
void SpectraSettings::set_x_label_offset(std::string x_label_offset) {
    this->x_label_offset = stod(x_label_offset);
}
void SpectraSettings::set_y_label_offset(std::string y_label_offset) {
    this->y_label_offset = stod(y_label_offset);
}


//--------------------------------------------------------------------------------------------------
// Default Constructor for PlotSettings
//--------------------------------------------------------------------------------------------------
PlotSettings::PlotSettings() {
    path_input_data = "output/output_trend.csv";
    path_output_figure = "output/output_trend.png";
    data_format = "xyy";
    title = "";
    x_label = "";
    y_label = "";
    y_min = 0; // if plotting detects max and min are the same, use default limits
    y_max = 0;
    x_min = 0; // if plotting detects max and min are the same, use default limits
    x_max = 0;
    x_res = 3200;
    y_res = 2400;
    x_log = 0;
    y_log = 0;
    x_grid = 0;
    y_grid = 0;
    x_num_divs = 0;
    y_num_divs = 0;
    legend_entries = {};
    // color_series = {"#000000","#C63822","#607FD5","#55A961"};
    //               black      blue       red      green     purple    grey      teal       salmon   light gr    fuchsia
    color_series = {"#000000","#4556d2","#C63822","#55A961","#75298e","#9ba5b1","#2fb5d4","#E79A9F","#4be0b0","#cb66ed"};
    // yellow to purple gradient:
    // color_series = {"#ffff00","#ece13c","#d9c452","#c5a661","#b18b6c","#9b6f75","#83547c","#693a81","#481f87","#00008b"};
    grayscale = 0;
    // color_error = {"#333333","#E79A9F","#6B8EF0","#69CF77"};
    // show_error;
    line_style = {1};
    line_width = {5};
    border_width = 5;
    legend = 1;
    legend_coords = {0.15,0.65,0.4,0.85};
    textbox = 0;
    textbox_coords = {0.15,0.4,0.4,0.6};
    textbox_text = {};
    plot_type = {"l"};
    legend_border_size = 0;
    legend_text_size = 0.035;
    marker_style = {20};
    marker_size = {5};
    margin_left = 0.1;
    margin_right = 0.1;
    margin_top = 0.1;
    margin_bottom = 0.1;
    x_label_offset = 1;
    y_label_offset = 0;
}

// Apply a value to a setting:
void PlotSettings::set_setting(std::string settings_name, std::string settings_value) {
    if (settings_value == ""){
        // If no setting value provided, do not apply anything
    }
    else if (settings_name == "path_input_data")
        this->set_path_input_data(settings_value);
    else if (settings_name == "path_output_figure")
        this->set_path_output_figure(settings_value);
    else if (settings_name == "data_format")
        this->set_data_format(settings_value);
    else if (settings_name == "title")
        this->set_title(settings_value);
    else if (settings_name == "x_label")
        this->set_x_label(settings_value);
    else if (settings_name == "y_label")
        this->set_y_label(settings_value);
    else if (settings_name == "x_min")
        this->set_x_min(settings_value);
    else if (settings_name == "x_max")
        this->set_x_max(settings_value);
    else if (settings_name == "y_min")
        this->set_y_min(settings_value);
    else if (settings_name == "y_max")
        this->set_y_max(settings_value);
    else if (settings_name == "x_res")
        this->set_x_res(settings_value);
    else if (settings_name == "y_res")
        this->set_y_res(settings_value);
    else if (settings_name == "x_log")
        this->set_x_log(settings_value);
    else if (settings_name == "y_log")
        this->set_y_log(settings_value);
    else if (settings_name == "x_grid")
        this->set_x_grid(settings_value);
    else if (settings_name == "y_grid")
        this->set_y_grid(settings_value);
    else if (settings_name == "y_num_divs")
        this->set_y_num_divs(settings_value);
    else if (settings_name == "x_num_divs")
        this->set_x_num_divs(settings_value);
    else if (settings_name == "legend_entries")
        this->set_legend_entries(settings_value);
    else if (settings_name == "color_series")
        this->set_color_series(settings_value);
    else if (settings_name == "grayscale")
        this->set_grayscale(settings_value);
    // else if (settings_name == "color_error")
    //     this->set_color_error(settings_value);
    // else if (settings_name == "show_error")
    //     this->set_show_error(settings_value);
    else if (settings_name == "line_style")
        this->set_line_style(settings_value);
    else if (settings_name == "line_width")
        this->set_line_width(settings_value);
    else if (settings_name == "border_width")
        this->set_border_width(settings_value);
    else if (settings_name == "legend")
        this->set_legend(settings_value);
    else if (settings_name == "legend_coords")
        this->set_legend_coords(settings_value);
    else if (settings_name == "textbox")
        this->set_textbox(settings_value);
    else if (settings_name == "textbox_coords")
        this->set_textbox_coords(settings_value);
    else if (settings_name == "textbox_text")
        this->set_textbox_text(settings_value);
    else if (settings_name == "plot_type")
        this->set_plot_type(settings_value);
    else if (settings_name == "legend_border_size")
        this->set_legend_border_size(settings_value);
    else if (settings_name == "legend_text_size")
        this->set_legend_text_size(settings_value);
    else if (settings_name == "marker_style")
        this->set_marker_style(settings_value);
    else if (settings_name == "marker_size")
        this->set_marker_size(settings_value);
    else if (settings_name == "margin_left")
        this->set_margin_left(settings_value);
    else if (settings_name == "margin_right")
        this->set_margin_right(settings_value);
    else if (settings_name == "margin_top")
        this->set_margin_top(settings_value);
    else if (settings_name == "margin_bottom")
        this->set_margin_bottom(settings_value);
    else if (settings_name == "x_label_offset")
        this->set_x_label_offset(settings_value);
    else if (settings_name == "y_label_offset")
        this->set_y_label_offset(settings_value);
    else
        throw std::logic_error("Unrecognized setting: " + settings_name 
            + ". Please refer to the README for allowed settings");
}

// Setter functions
void PlotSettings::set_path_input_data(std::string path_input_data) {
    this->path_input_data = path_input_data;
}
void PlotSettings::set_path_output_figure(std::string path_output_figure) {
    this->path_output_figure = path_output_figure;
}
void PlotSettings::set_data_format(std::string data_format) {
    this->data_format = data_format;
}
void PlotSettings::set_title(std::string title) {
    this->title = title;
}
void PlotSettings::set_x_label(std::string x_label) {
    this->x_label = x_label;
}
void PlotSettings::set_y_label(std::string y_label) {
    this->y_label = y_label;
}
void PlotSettings::set_x_min(std::string x_min) {
    this->x_min = stod(x_min);
}
void PlotSettings::set_x_max(std::string x_max) {
    this->x_max = stod(x_max);
}
void PlotSettings::set_y_min(std::string y_min) {
    this->y_min = stod(y_min);
}
void PlotSettings::set_y_max(std::string y_max) {
    this->y_max = stod(y_max);
}
void PlotSettings::set_x_res(std::string x_res) {
    this->x_res = stoi(x_res);
}
void PlotSettings::set_y_res(std::string y_res) {
    this->y_res = stoi(y_res);
}
void PlotSettings::set_x_log(std::string x_log) {
    this->x_log = stoi(x_log);
}
void PlotSettings::set_y_log(std::string y_log) {
    this->y_log = stoi(y_log);
}
void PlotSettings::set_x_grid(std::string x_grid) {
    this->x_grid = stoi(x_grid);
}
void PlotSettings::set_y_grid(std::string y_grid) {
    this->y_grid = stoi(y_grid);
}
void PlotSettings::set_x_num_divs(std::string x_num_divs) {
    this->x_num_divs = stoi(x_num_divs);
}
void PlotSettings::set_y_num_divs(std::string y_num_divs) {
    this->y_num_divs = stoi(y_num_divs);
}
void PlotSettings::set_legend_entries(std::string legend_entries) {
    stringToSVector(legend_entries,this->legend_entries);
}
void PlotSettings::set_color_series(std::string color_series) {
    stringToSVector(color_series,this->color_series);
}
void PlotSettings::set_grayscale(std::string grayscale) {
    this->grayscale = stoi(grayscale);
} 
// void PlotSettings::set_color_error(std::string color_error) {
//     stringToSVector(color_error,this->color_error);
// }
// void PlotSettings::set_show_error(std::string show_error) {
//     stringToIVector(show_error,this->show_error);
// }
void PlotSettings::set_line_style(std::string line_style) {
    stringToIVector(line_style,this->line_style);
}
void PlotSettings::set_line_width(std::string line_width) {
    stringToIVector(line_width,this->line_width);
}
void PlotSettings::set_border_width(std::string border_width) {
    this->border_width = stoi(border_width);
}
void PlotSettings::set_legend(std::string legend) {
    this->legend = stoi(legend);
} 
void PlotSettings::set_legend_coords(std::string legend_coords) {
    stringToDVector(legend_coords,this->legend_coords);
}
void PlotSettings::set_textbox(std::string textbox) {
    this->textbox = stoi(textbox);
} 
void PlotSettings::set_textbox_coords(std::string textbox_coords) {
    stringToDVector(textbox_coords,this->textbox_coords);
}
void PlotSettings::set_textbox_text(std::string textbox_text) {
    stringToSVector(textbox_text,this->textbox_text);
}
void PlotSettings::set_plot_type(std::string plot_type) {
    stringToSVector(plot_type,this->plot_type);
}
void PlotSettings::set_legend_border_size(std::string legend_border_size) {
    this->legend_border_size = stoi(legend_border_size);
}
void PlotSettings::set_legend_text_size(std::string legend_text_size) {
    this->legend_text_size = stod(legend_text_size);
}
void PlotSettings::set_marker_style(std::string marker_style) {
    stringToIVector(marker_style,this->marker_style);
}
void PlotSettings::set_marker_size(std::string marker_size) {
    stringToIVector(marker_size,this->marker_size);
}
void PlotSettings::set_margin_left(std::string margin_left) {
    this->margin_left = stod(margin_left);
}
void PlotSettings::set_margin_right(std::string margin_right) {
    this->margin_right = stod(margin_right);
}
void PlotSettings::set_margin_top(std::string margin_top) {
    this->margin_top = stod(margin_top);
}
void PlotSettings::set_margin_bottom(std::string margin_bottom) {
    this->margin_bottom = stod(margin_bottom);
}
void PlotSettings::set_x_label_offset(std::string x_label_offset) {
    this->x_label_offset = stod(x_label_offset);
}
void PlotSettings::set_y_label_offset(std::string y_label_offset) {
    this->y_label_offset = stod(y_label_offset);
}


//--------------------------------------------------------------------------------------------------
// Default Constructor for SurfaceSettings
//--------------------------------------------------------------------------------------------------
SurfaceSettings::SurfaceSettings() {
    path_input_data = "output/output_surface.csv";
    path_output_figure = "output/output_sruface.png";
    title = "";
    x_label = "";
    y_label = "";
    z_label = "";
    y_min = 0; // if plotting detects max and min are the same, use default limits
    y_max = 0;
    x_min = 0; // if plotting detects max and min are the same, use default limits
    x_max = 0;
    z_min = 0; // if plotting detects max and min are the same, use default limits
    z_max = 0;
    x_res = 800;
    y_res = 600;
    x_num_divs = 206;
    z_num_divs = 203;
    color_palette = 55;
    num_color_bins = 40;
    border_width = 1;
}

// Apply a value to a setting:
void SurfaceSettings::set_setting(std::string settings_name, std::string settings_value) {
    if (settings_value == ""){
        // If no setting value provided, do not apply anything
    }
    else if (settings_name == "path_input_data")
        this->set_path_input_data(settings_value);
    else if (settings_name == "path_output_figure")
        this->set_path_output_figure(settings_value);
    else if (settings_name == "title")
        this->set_title(settings_value);
    else if (settings_name == "x_label")
        this->set_x_label(settings_value);
    else if (settings_name == "y_label")
        this->set_y_label(settings_value);
    else if (settings_name == "z_label")
        this->set_z_label(settings_value);
    else if (settings_name == "x_min")
        this->set_x_min(settings_value);
    else if (settings_name == "x_max")
        this->set_x_max(settings_value);
    else if (settings_name == "y_min")
        this->set_y_min(settings_value);
    else if (settings_name == "y_max")
        this->set_y_max(settings_value);
    else if (settings_name == "z_min")
        this->set_z_min(settings_value);
    else if (settings_name == "z_max")
        this->set_z_max(settings_value);
    else if (settings_name == "x_res")
        this->set_x_res(settings_value);
    else if (settings_name == "y_res")
        this->set_y_res(settings_value);
    else if (settings_name == "z_num_divs")
        this->set_z_num_divs(settings_value);
    else if (settings_name == "x_num_divs")
        this->set_x_num_divs(settings_value);
    else if (settings_name == "color_palette")
        this->set_color_palette(settings_value);
    else if (settings_name == "num_color_bins")
        this->set_num_color_bins(settings_value);
    else if (settings_name == "border_width")
        this->set_border_width(settings_value);

    else
        throw std::logic_error("Unrecognized setting: " + settings_name 
            + ". Please refer to the README for allowed settings");
}

// Setter functions
void SurfaceSettings::set_path_input_data(std::string path_input_data) {
    this->path_input_data = path_input_data;
}
void SurfaceSettings::set_path_output_figure(std::string path_output_figure) {
    this->path_output_figure = path_output_figure;
}
void SurfaceSettings::set_title(std::string title) {
    this->title = title;
}
void SurfaceSettings::set_x_label(std::string x_label) {
    this->x_label = x_label;
}
void SurfaceSettings::set_y_label(std::string y_label) {
    this->y_label = y_label;
}
void SurfaceSettings::set_z_label(std::string z_label) {
    this->z_label = z_label;
}
void SurfaceSettings::set_x_min(std::string x_min) {
    this->x_min = stoi(x_min);
}
void SurfaceSettings::set_x_max(std::string x_max) {
    this->x_max = stoi(x_max);
}
void SurfaceSettings::set_y_min(std::string y_min) {
    this->y_min = stod(y_min);
}
void SurfaceSettings::set_y_max(std::string y_max) {
    this->y_max = stod(y_max);
}
void SurfaceSettings::set_z_min(std::string z_min) {
    this->z_min = stod(z_min);
}
void SurfaceSettings::set_z_max(std::string z_max) {
    this->z_max = stod(z_max);
}
void SurfaceSettings::set_x_res(std::string x_res) {
    this->x_res = stoi(x_res);
}
void SurfaceSettings::set_y_res(std::string y_res) {
    this->y_res = stoi(y_res);
}
void SurfaceSettings::set_x_num_divs(std::string x_num_divs) {
    this->x_num_divs = stoi(x_num_divs);
}
void SurfaceSettings::set_z_num_divs(std::string z_num_divs) {
    this->z_num_divs = stoi(z_num_divs);
}
void SurfaceSettings::set_color_palette(std::string color_palette) {
    this->color_palette = stoi(color_palette);
}
void SurfaceSettings::set_num_color_bins(std::string num_color_bins) {
    this->num_color_bins = stoi(num_color_bins);
}
void SurfaceSettings::set_border_width(std::string border_width) {
    this->border_width = stoi(border_width);
}