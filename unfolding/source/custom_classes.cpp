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
    norm = 1.21;
    error = 0.0;
    f_factor = 7.0;
    cutoff = 10000;
    uncertainty_type = "poisson";
    num_poisson_samples = 50;
    meas_units = "nc";
    // MAP specific
    beta = 0.0;
    prior = "mrp";
    // MLEM-STOP specific
    cps_crossover = 40000;
    sigma_j=0.5;
    // Optimize specific
    min_num_iterations = 500;
    max_num_iterations = 10000;
    iteration_increment = 50;
    min_beta = 0.0;
    max_beta = 1.0;
    parameter_of_interest = "fluence";
    algorithm = "mlem";
    trend_type = "ratio";
    path_output_spectra = "output/output_spectra.csv";
    generate_report = 1;
    path_report = "";
    generate_figure = 0;
    path_figure = "";
    auto_output_path = "output/auto.csv";
    derivatives = 0;
    measurements_path = "input/measurements.txt";
    input_spectrum_path = "input/spectrum_step.csv";
    energy_bins_path = "input/energy_bins.csv";
    system_response_path = "input/nns_response.csv";
    icrp_factors_path = "input/icrp_conversions.csv";
    ref_spectrum_path = "";
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
    else if (settings_name == "num_poisson_samples")
        this->set_num_poisson_samples(atoi(settings_value.c_str()));
    else if (settings_name == "meas_units")
        this->set_meas_units(settings_value);
    else if (settings_name == "beta")
        this->set_beta(atof(settings_value.c_str()));
    else if (settings_name == "prior")
        this->set_prior(settings_value);
    else if (settings_name == "cps_crossover")
        this->set_cps_crossover(atoi(settings_value.c_str()));
    else if (settings_name == "sigma_j")
        this->set_sigma_j(atof(settings_value.c_str()));
    else if (settings_name == "min_num_iterations")
        this->set_min_num_iterations(atoi(settings_value.c_str()));
    else if (settings_name == "max_num_iterations")
        this->set_max_num_iterations(atoi(settings_value.c_str()));
    else if (settings_name == "iteration_increment")
        this->set_iteration_increment(atoi(settings_value.c_str()));
    else if (settings_name == "min_beta")
        this->set_min_beta(atof(settings_value.c_str()));
    else if (settings_name == "max_beta")
        this->set_max_beta(atof(settings_value.c_str()));
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
    else if (settings_name == "auto_output_path")
        this->set_auto_output_path(settings_value);
    else if (settings_name == "derivatives")
        this->set_derivatives(atoi(settings_value.c_str()));
    else if (settings_name == "measurements_path")
        this->set_measurements_path(settings_value);
    else if (settings_name == "input_spectrum_path")
        this->set_input_spectrum_path(settings_value);
    else if (settings_name == "energy_bins_path")
        this->set_energy_bins_path(settings_value);
    else if (settings_name == "system_response_path")
        this->set_system_response_path(settings_value);
    else if (settings_name == "icrp_factors_path")
        this->set_icrp_factors_path(settings_value);
    else if (settings_name == "ref_spectrum_path")
        this->set_ref_spectrum_path(settings_value);
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
void UnfoldingSettings::set_num_poisson_samples(int num_poisson_samples) {
    this->num_poisson_samples = num_poisson_samples;
}
void UnfoldingSettings::set_meas_units(std::string meas_units) {
    this->meas_units = meas_units;
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
void UnfoldingSettings::set_auto_output_path(std::string auto_output_path) {
    this->auto_output_path = auto_output_path;
}
void UnfoldingSettings::set_derivatives(int derivatives) {
    this->derivatives = derivatives;
}
void UnfoldingSettings::set_measurements_path(std::string measurements_path) {
    this->measurements_path = measurements_path;
}
void UnfoldingSettings::set_input_spectrum_path(std::string input_spectrum_path) {
    this->input_spectrum_path = input_spectrum_path;
}
void UnfoldingSettings::set_energy_bins_path(std::string energy_bins_path) {
    this->energy_bins_path = energy_bins_path;
}
void UnfoldingSettings::set_system_response_path(std::string system_response_path) {
    this->system_response_path = system_response_path;
}
void UnfoldingSettings::set_icrp_factors_path(std::string icrp_factors_path) {
    this->icrp_factors_path = icrp_factors_path;
}
void UnfoldingSettings::set_ref_spectrum_path(std::string ref_spectrum_path) {
    this->ref_spectrum_path = ref_spectrum_path;
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
void UnfoldingReport::set_num_poisson_samples(int num_poisson_samples) {
    this->num_poisson_samples = num_poisson_samples;
}
void UnfoldingReport::set_git_commit(std::string git_commit) {
    this->git_commit = git_commit;
}
void UnfoldingReport::set_measurements(std::vector<double>& measurements) {
    this->measurements = measurements;
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
    rfile << std::left << std::setw(sw) << "Git commit number: " << git_commit << "\n";
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
    rfile << std::left << std::setw(sw) << "Number of poisson samples:" << num_poisson_samples << "\n";
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
    if (meas_units == "cps")
        units_string = "CPS";
    else
        units_string = "Charge (nC)";
    
    rfile << "Measurement\n\n";
    rfile << std::left << std::setw(sw) << "Delivered dose:" << dose_mu << " MU\n";
    rfile << std::left << std::setw(sw) << "Delivered doserate:" << doserate_mu << " MU/min\n";
    rfile << std::left << std::setw(sw) << "Measurement duration:" << duration << " s\n\n";
    // rfile << "Measured Data (measurement duration: " << duration << "s)\n\n";
    rfile << std::left << std::setw(cw) << "# of moderators" << units_string << "\n";
    rfile << std::left << std::setw(cw) << COLSTRING << COLSTRING << "\n";
    for (int i=0; i<num_measurements; i++) {
        rfile << std::left << std::setw(cw) << num_measurements-1-i << measurements[i] << "\n";
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
    if (algorithm == "mlemstop") {
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
    input_filename = "output_spectra.csv";
    input_dir = "output/";
    output_filename = "neutron_spectra.png";
    output_dir = "output/";
    title = "";
    x_label = "Energy (MeV)";
    y_label = "Fluence (n #upoint cm^{-2} s^{-1})";
    y_min = 0; // if plotting detects max and min are the same, use default limits
    y_max = 0;
    x_res = 800;
    y_res = 600;
    y_num_divs = 0;
    legend_entries = {};
    color_series = {"#000000","#C63822","#607FD5","#55A961","#75298e","#ce884a","#33889b"};
    color_error = {"#333333","#E79A9F","#6B8EF0","#69CF77","#ad4ace","#e2a56f","#69c6db"};
    grayscale = 0;
    show_error = {1};
    error_style = "E2";
    rows_per_spectrum = 3;
    line_style = {1};
    line_width = {1};
    border_width = 1;
    legend = 1;
    legend_coords = {0.15,0.75,0.4,0.85};
    textbox = 0;
    textbox_coords = {0.15,0.4,0.4,0.6};
    textbox_text = {};
    plot_per_mu = 0;
    number_mu = {0};
    duration = {0};
    normalize = 0;
}

// Apply a value to a setting:
void SpectraSettings::set_setting(std::string settings_name, std::string settings_value) {
    if (settings_value == ""){
        // If no setting value provided, do not apply anything
    }
    else if (settings_name == "input_filename")
        this->set_input_filename(settings_value);
    else if (settings_name == "input_dir")
        this->set_input_dir(settings_value);
    else if (settings_name == "output_filename")
        this->set_output_filename(settings_value);
    else if (settings_name == "output_dir")
        this->set_output_dir(settings_value);
    else if (settings_name == "title")
        this->set_title(settings_value);
    else if (settings_name == "x_label")
        this->set_x_label(settings_value);
    else if (settings_name == "y_label")
        this->set_y_label(settings_value);
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
    // else
    //     throw std::logic_error("Unrecognized setting: " + settings_name 
    //      + ". Please refer to the README for allowed settings");
}

void SpectraSettings::set_input_filename(std::string input_filename) {
    this->input_filename = input_filename;
}
void SpectraSettings::set_input_dir(std::string input_dir) {
    this->input_dir = input_dir;
}
void SpectraSettings::set_output_filename(std::string output_filename) {
    this->output_filename = output_filename;
}
void SpectraSettings::set_output_dir(std::string output_dir) {
    this->output_dir = output_dir;
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


//--------------------------------------------------------------------------------------------------
// Default Constructor for PlotSettings
//--------------------------------------------------------------------------------------------------
PlotSettings::PlotSettings() {
    input_filename = "poi_output_mlem.csv";
    input_dir = "output/";
    output_filename = "poi_output_mlem.png";
    output_dir = "output/";
    data_format = "xyy";
    title = "";
    x_label = "";
    y_label = "";
    y_min = 0; // if plotting detects max and min are the same, use default limits
    y_max = 0;
    x_min = 0; // if plotting detects max and min are the same, use default limits
    x_max = 0;
    x_res = 800;
    y_res = 600;
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
    line_width = {2};
    border_width = 1;
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
    else if (settings_name == "input_filename")
        this->set_input_filename(settings_value);
    else if (settings_name == "input_dir")
        this->set_input_dir(settings_value);
    else if (settings_name == "output_filename")
        this->set_output_filename(settings_value);
    else if (settings_name == "output_dir")
        this->set_output_dir(settings_value);
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
void PlotSettings::set_input_filename(std::string input_filename) {
    this->input_filename = input_filename;
}
void PlotSettings::set_input_dir(std::string input_dir) {
    this->input_dir = input_dir;
}
void PlotSettings::set_output_filename(std::string output_filename) {
    this->output_filename = output_filename;
}
void PlotSettings::set_output_dir(std::string output_dir) {
    this->output_dir = output_dir;
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
    input_filename = "poi_output_map.csv";
    input_dir = "output/";
    output_filename = "poi_output_map.png";
    output_dir = "output/";
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
    else if (settings_name == "input_filename")
        this->set_input_filename(settings_value);
    else if (settings_name == "input_dir")
        this->set_input_dir(settings_value);
    else if (settings_name == "output_filename")
        this->set_output_filename(settings_value);
    else if (settings_name == "output_dir")
        this->set_output_dir(settings_value);
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
void SurfaceSettings::set_input_filename(std::string input_filename) {
    this->input_filename = input_filename;
}
void SurfaceSettings::set_input_dir(std::string input_dir) {
    this->input_dir = input_dir;
}
void SurfaceSettings::set_output_filename(std::string output_filename) {
    this->output_filename = output_filename;
}
void SurfaceSettings::set_output_dir(std::string output_dir) {
    this->output_dir = output_dir;
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