//**************************************************************************************************
// This program unfolds measurements (nA or cps) obtained using a Nested Neutron Spectrometer
// to calculate one of various parameters of interest (e.g. dose, fluence, MLEM ratio) with
// increasing iteration number. In the case of MAP unfolding, the the parameter is also calculated
// as a function of beta.
//
// Output:
//  - For MLEM a CSV file with the iteration numbers and corresponding POI values are output. If
//    the output file already exists, the current data will be appended.
//  - For MAP, a CSV file with the iteration numbers and corresponding POI values corresponding to
//    a range of beta values are ouput. Each row in the output file correspond to a particular beta
//    value. 
//**************************************************************************************************

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdlib.h>
#include <vector>

// Local
#include "custom_classes.h"
#include "fileio.h"
#include "handle_args.h"
#include "root_helpers.h"
#include "physics_calculations.h"

int main(int argc, char* argv[])
{
    // Put arguments in vector for easier processing
    std::vector<std::string> arg_vector;
    for (int i = 1; i < argc; i++) {
        arg_vector.push_back(argv[i]);
    }

    // NOTE: Indices are linked between the following arrays and vectors (i.e. input_files[0]
    // corresponds to input_file_flags[0] and input_file_defaults[0])
    // Array that stores the allowed options that specify input files
    // Add new options at end of array
    const int num_ifiles = 1;
    std::string input_file_flags_arr[num_ifiles] = {
        "--configuration"
    };
    // Array that stores default filename for each input file
    std::string input_file_defaults_arr[num_ifiles] = {
        "input/unfold_trend.cfg"
    };

    // Convert arrays to vectors b/c easier to work with
    std::vector<std::string> input_files; // Store the actual input filenames to be used
    std::vector<std::string> input_file_flags;
    std::vector<std::string> input_file_defaults;
    for (int i=0; i<num_ifiles; i++) {
        input_files.push_back("");
        input_file_flags.push_back(input_file_flags_arr[i]);
        input_file_defaults.push_back(input_file_defaults_arr[i]);
    }

    // Use provided arguments (files) and/or defaults to determine the input files to be used
    for (int i=0; i<num_ifiles; i++) {
        setfile(arg_vector, input_file_flags[i], input_file_defaults[i], input_files[i]);
    }

    // Notify user if unknown parameters were received
    checkUnknownParameters(arg_vector, input_file_flags);

    // Apply some settings read in from a config file
    UnfoldingSettings settings;
    setSettings(input_files[0], settings);

    settings.set_f_factor(settings.f_factor / 1e6); // Convert f_factor from fA/cps to nA/cps

    // Read in measurements from file
    std::vector<double> measurements_nc;
    std::vector<double> measurements;
    int num_measurements = 0;
    measurements = getMeasurements(settings);
    num_measurements = measurements.size();
    std::reverse(measurements.begin(),measurements.end()); // readin 7-0 but want 0-7

    // Handle nC input (save nC values for report, then convert to CPS)
    if (settings.meas_units == "nc") {
        measurements_nc = measurements;
        for (int i_meas=0; i_meas < num_measurements; i_meas++) {
            measurements[i_meas] = measurements[i_meas]*settings.norm/settings.f_factor/settings.duration;
        }
    }

    // Process data wherein multiple measurements were acquired for each shell
    std::vector<double> std_errors;
    if (settings.num_meas_per_shell > 1) {
        processMeasurements(num_measurements,settings.num_meas_per_shell,measurements,std_errors);
        num_measurements = num_measurements / settings.num_meas_per_shell;
    }
    // Do not allow use of gaussian sampling technique if only one measurement per shell is
    // provided, as the standard deviation is unknown.
    else if (settings.num_meas_per_shell == 1 && settings.uncertainty_type == "gaussian"){
        throw std::logic_error("Cannot generate Gaussian-sampled pseudo-measurements with only single measurement per shell.");
    }
    // Require at least 1 measurement per shell
    else if (settings.num_meas_per_shell < 1) { // if settings.num_meas_per_shell is 0 or negative
        throw std::logic_error("Number of measurements per shell must be >= 1");
    }

    // std::vector<double> measurements_nc = getMeasurements(input_files[0], irradiation_conditions, 
    //    dose_mu, doserate_mu, duration);
    // int num_measurements = measurements_nc.size();

    // // Convert measured charge in nC to counts per second
    // // Re-order measurments from (7 moderators to 0) to (0 moderators to 7)
    // std::vector<double> measurements;
    // for (int index=0; index < num_measurements; index++) {
    //     double measurement_cps = measurements_nc[num_measurements-index-1]*settings.norm/settings.f_factor/duration;
    //     measurements.push_back(measurement_cps);
    // }

    //----------------------------------------------------------------------------------------------
    // Generate the energy bins matrix:
    //  - size = # of energy bins
    // Input the energies from energy bins file:
    //  - values in units of [MeV]
    //----------------------------------------------------------------------------------------------
    std::vector<double> energy_bins;
    readInputFile1D(settings.path_energy_bins,energy_bins);

    int num_bins = energy_bins.size();

    //----------------------------------------------------------------------------------------------
    // Generate the detector response matrix (representing the dector response function):
    //  - outer size = # of measurements
    //  - inner size = # of energy bins
    // Detector response value for each # of moderators for each energy (currently 52 energies)
    // Input the response functions
    //  - values in units of [cm^2]
    //
    // The response function accounts for variable number of (n,p) reactions in He-3 for each
    // moderators, as a function of energy. Calculated by vendor using MC
    //----------------------------------------------------------------------------------------------
    std::vector<std::vector<double>> nns_response;
    readInputFile2D(settings.path_system_response,nns_response);
    checkDimensions(num_measurements, "number of measurements", nns_response.size(), "NNS response");
    checkDimensions(num_bins, "number of energy bins", nns_response[0].size(), "NNS response");

    //----------------------------------------------------------------------------------------------
    // Generate the inital spectrum matrix to input into the unfolding algorithm:
    //  - size = # of energy bins
    // Input from the input spectrum file
    //  - values are neutron fluence rates [neutrons cm^-2 s^-1])
    //  - Currently (2017-08-16) input a step function (high at thermals & lower), because a flat 
    //  spectrum underestimates (does not yield any) thermal neutrons
    //----------------------------------------------------------------------------------------------
    std::vector<double> initial_spectrum;
    readInputFile1D(settings.path_input_spectrum,initial_spectrum);
    checkDimensions(num_bins, "number of energy bins", initial_spectrum.size(), "Input spectrum");

    std::vector<double> spectrum = initial_spectrum; // save the initial spectrum for report output

    //----------------------------------------------------------------------------------------------
    // Generate the ICRP conversion matrix (factors to convert fluence to ambient dose equivalent):
    //  - size = # of energy bins
    // Input from icrp factors file:
    //  - values are in units of [pSv cm^2]
    //  - H values were obtained by linearly interopolating tabulated data to match energy bins used
    // Page 200 of document (ICRP 74 - ATables.pdf)
    //----------------------------------------------------------------------------------------------
    std::vector<double> icrp_factors;
    readInputFile1D(settings.path_icrp_factors,icrp_factors);
    checkDimensions(num_bins, "number of energy bins", icrp_factors.size(), "Number of ICRP factors");

    //----------------------------------------------------------------------------------------------
    // Run the automatic unfolding algorithm.
    // Note: the normalized system matrix is calculated first. It is required in unfolding, and is
    // a constant value.
    //----------------------------------------------------------------------------------------------
    std::vector<double> normalized_response = normalizeResponse(num_bins, num_measurements, nns_response);

    std::vector<double> mlem_ratio; // vector that stores the ratio between measured data and MLEM estimated data
    std::vector<double> mlem_correction; // vector that stores the ratio between measured data and MLEM estimated data
    std::vector<double> mlem_estimate;

    //----------------------------------------------------------------------------------------------
    // Output correction factors (52 values applied to spectrum, NOT to measurements).
    // Each row of output file has the 52 correction factors for a specific N
    // Visualize with plot_lines (logarithmic x-axis)
    //----------------------------------------------------------------------------------------------
    if (settings.algorithm == "correction_factors") {
        // Create vector of number of iterations
        int num_increments = ((settings.iteration_max - settings.iteration_min) / settings.iteration_increment)+1;
        std::vector<int> num_iterations_vector = linearSpacedIntegerVector(
            settings.iteration_min,settings.iteration_max,num_increments
        );
        int num_iteration_samples = num_iterations_vector.size();

        std::vector<double> current_spectrum = initial_spectrum; // the reconstructed spectrum

        // Add the energy bins to the file
        std::ostringstream results_stream;
        results_stream << "Energy (MeV),";
        for (int i_bin = 0; i_bin < num_bins; i_bin++) {
            results_stream << energy_bins[i_bin];
            if (i_bin != num_bins-1)
                results_stream << ",";
        }
        results_stream << "\n";

        int total_num_iterations = 0;
        for (int i_num=0; i_num < num_iteration_samples; i_num++) {
            int num_iterations;

            // only do # of iterations since previous run (i.e. don't do 3000, then 4000. Do
            // 3000 then iteration_size more, etc.)
            if (i_num == 0)
                num_iterations = num_iterations_vector[i_num];
            else
                num_iterations = num_iterations_vector[i_num]-num_iterations_vector[i_num-1];
            runMLEM(num_iterations, settings.error, num_measurements, num_bins, measurements, 
                current_spectrum, nns_response, normalized_response, mlem_ratio, mlem_correction, mlem_estimate
            );

            total_num_iterations += num_iterations;

            // Add the reconstructed measured data for current number of iterations to the file
            // Add the correction factors to the file
            results_stream << "k = " << total_num_iterations << ",";
            for (int i_bin = 0; i_bin < num_bins; i_bin++) {
                results_stream << mlem_correction[i_bin];

                if (i_bin != num_bins-1)
                    results_stream << ",";
            }
            results_stream << "\n";
        }        

        // Save results for parameter of interest to CSV file
        std::ofstream output_file;
        output_file.open(settings.path_output_trend, std::ios_base::out);
        std::string results_string = results_stream.str();
        output_file << results_string;
        output_file.close();

        std::cout << "Saved correction factors to " << settings.path_output_trend << "\n";
    }

    //----------------------------------------------------------------------------------------------
    // Output reconstructed measured data (8 values) either as absolute values (CPS) or as ratio
    // with real measured data.
    // First row contains number of moderators
    // First row contains the real measured data (either cps or ratio=1)
    // Remaining rows have the reconstructed measurements (either cps or ratio) for a specific N
    // Visualize with plot_lines (scatter plots)
    //----------------------------------------------------------------------------------------------
    if (settings.algorithm == "trend") {
        // Create vector of number of iterations
        int num_increments = ((settings.iteration_max - settings.iteration_min) / settings.iteration_increment)+1;
        std::vector<int> num_iterations_vector = linearSpacedIntegerVector(
            settings.iteration_min,settings.iteration_max,num_increments
        );
        int num_iteration_samples = num_iterations_vector.size();

        std::vector<double> current_spectrum = initial_spectrum; // the reconstructed spectrum

        // Add the Moderator numbers to the file
        std::ostringstream results_stream;
        results_stream << "Number of moderators,";
        for (int i_meas = 0; i_meas < num_measurements; i_meas++) {
            results_stream << i_meas;
            if (i_meas != num_measurements-1)
                results_stream << ",";
        }
        results_stream << "\n";

        // Add the measured data to the file
        results_stream << "Measured data,";
        for (int i_meas = 0; i_meas < num_measurements; i_meas++) {
            // If comparing CPS data:
            if (settings.trend_type == "cps") {
                results_stream << measurements[i_meas];
            }
            // If comparing the reconstructed MLEM ratios:
            else if (settings.trend_type == "ratio") {
                results_stream << 1; // ratio of measured data and itself is 1
            }


            if (i_meas != num_measurements-1)
                results_stream << ",";
        }
        results_stream << "\n";

        int total_num_iterations = 0;
        for (int i_num=0; i_num < num_iteration_samples; i_num++) {
            int num_iterations;

            // only do # of iterations since previous run (i.e. don't do 3000, then 4000. Do
            // 3000 then iteration_size more, etc.)
            if (i_num == 0)
                num_iterations = num_iterations_vector[i_num];
            else
                num_iterations = num_iterations_vector[i_num]-num_iterations_vector[i_num-1];
            runMLEM(num_iterations, settings.error, num_measurements, num_bins, measurements, 
                current_spectrum, nns_response, normalized_response, mlem_ratio, mlem_correction, mlem_estimate
            );

            total_num_iterations += num_iterations;
            // Add the reconstructed measured data for current number of iterations to the file
            results_stream << "N = " << total_num_iterations << ",";
            for (int i_meas = 0; i_meas < num_measurements; i_meas++) {
                // If comparing CPS data, use the returned ratio to calculate the reconstructed measurement:
                if (settings.trend_type == "cps") {
                    double recon_meas = measurements[i_meas] / mlem_ratio[i_meas];
                    results_stream << recon_meas;
                }
                // If comparing the reconstructed MLEM ratios:
                else if (settings.trend_type == "ratio") {
                    results_stream << mlem_ratio[i_meas];
                }


                if (i_meas != num_measurements-1)
                    results_stream << ",";
            }
            results_stream << "\n";
        }        

        // Save results for parameter of interest to CSV file
        std::ofstream output_file;
        output_file.open(settings.path_output_trend, std::ios_base::out);
        std::string results_string = results_stream.str();
        output_file << results_string;
        output_file.close();

        std::cout << "Saved reconstruced measured data to " << settings.path_output_trend << "\n";
    }

    //----------------------------------------------------------------------------------------------
    // Calculate some parameter of interest (POI) at specified numbers of iterations (N). Iterate
    // through N values as indicated by user. Append results to existing file if it exists. 
    // Output:
    // First row contains the iteration numbers (only done if output file is empty)
    // Other row contains the POI value at each N. A single execution of this program will only add
    //  one result line to the the output file. But can run mulitple times for different measured
    //  data to append to existing file
    // Visualize with plot_lines
    //----------------------------------------------------------------------------------------------
    if (settings.algorithm == "mlem") {
        // Create vector of number of iterations
        int num_increments = ((settings.iteration_max - settings.iteration_min) / settings.iteration_increment)+1;
        std::vector<int> num_iterations_vector = linearSpacedIntegerVector(
            settings.iteration_min,settings.iteration_max,num_increments
        );
        int num_iteration_samples = num_iterations_vector.size();

        std::vector<double> current_spectrum = initial_spectrum; // the reconstructed spectrum

        std::vector<double> derivative_poi_values; // hold poi_values in vector if doing derivative calculations

        // Create stream to append results. First row is number of iteration increments
        // determine if file exists
        // std::string map_filename = "output/poi_output_mlem.csv";
        std::ifstream rfile(settings.path_output_trend);
        bool file_empty = is_empty(rfile);
        // bool file_exists = rfile.good();
        rfile.close();

        // If the file is empty, make the first line the number of iterations
        std::ostringstream results_stream;
        if (file_empty) {
            results_stream << "Number of iterations,";
            for (int i_num = 0; i_num < num_iteration_samples; i_num++) {
                results_stream << num_iterations_vector[i_num];
                if (i_num != num_iteration_samples-1)
                    results_stream << ",";
            }
            results_stream << "\n";
        }
        // Used if calculating RMSD
        std::vector<double> ref_spectrum;
        if (settings.parameter_of_interest == "rms" || settings.parameter_of_interest == "nrmsd" 
            || settings.parameter_of_interest == "chi_squared_g") 
        {
            readInputFile1D(settings.path_ref_spectrum,ref_spectrum);
        }

        // Loop through number of iterations
        for (int i_num=0; i_num < num_iteration_samples; i_num++) {
            int num_iterations;

            // only do # of iterations since previous run (i.e. don't do 3000, then 4000. Do
            // 3000 then iteration_size more, etc.)
            if (i_num == 0)
                num_iterations = num_iterations_vector[i_num];
            else
                num_iterations = num_iterations_vector[i_num]-num_iterations_vector[i_num-1];
            runMLEM(num_iterations, settings.error, num_measurements, num_bins, measurements, current_spectrum, 
                nns_response, normalized_response, mlem_ratio, mlem_correction, mlem_estimate
            );
        
            // Calculate one of the following parameters of interest & save to results stream
            double poi_value = 0;
            // Total fluence
            if (settings.parameter_of_interest == "total_fluence") {
                if (i_num == 0)
                    results_stream << settings.irradiation_conditions << ",";
                poi_value = calculateTotalFlux(num_bins,current_spectrum);
            }
            // Total dose:
            else if (settings.parameter_of_interest == "total_dose") {
                if (i_num == 0)
                    results_stream << "Total dose,";
                poi_value = calculateDose(num_bins, current_spectrum, icrp_factors);
            }
            // Maximum "error" in MLEM ratio
            else if (settings.parameter_of_interest == "max_mlem_ratio") {
                if (i_num == 0)
                    results_stream << settings.irradiation_conditions << ",";
                poi_value = calculateMaxRatio(num_measurements,mlem_ratio);
            }
            // Average "error" in MLEM ratio
            else if (settings.parameter_of_interest == "avg_mlem_ratio") {
                if (i_num == 0)
                    results_stream << settings.irradiation_conditions << ",";
                poi_value = calculateAvgRatio(num_measurements,mlem_ratio);
            }
            else if (settings.parameter_of_interest == "j_factor") {
                if (i_num == 0)
                    results_stream << settings.irradiation_conditions << ",";
                poi_value = calculateJFactor(num_measurements,measurements,mlem_estimate);
            }
            else if (settings.parameter_of_interest == "j_factor2") {
                if (i_num == 0)
                    results_stream << settings.irradiation_conditions << ",";
                poi_value = calculateJFactor2(num_measurements,measurements,mlem_estimate);
            }
            else if (settings.parameter_of_interest == "noise") {
                if (i_num == 0)
                    results_stream << settings.irradiation_conditions << ",";
                poi_value = calculateNoise(15,30,current_spectrum);
            }
            else if (settings.parameter_of_interest == "reduced_chi_squared") {
                if (i_num == 0)
                    results_stream << settings.irradiation_conditions << ",";
                poi_value = calculateChiSquared(i_num,num_bins,num_measurements,current_spectrum,measurements,mlem_ratio);
            }
            else if (settings.parameter_of_interest == "rms") {
                if (i_num == 0)
                    results_stream << settings.irradiation_conditions << ",";
                // std::vector<double> normalized_spectrum = normalizeVector(current_spectrum);
                std::vector<double> normalized_spectrum = current_spectrum;
                poi_value = calculateRMSEstimator(num_bins,ref_spectrum,normalized_spectrum);
            }
            else if (settings.parameter_of_interest == "nrmsd") {
                if (i_num == 0)
                    results_stream << settings.irradiation_conditions << ",";
                // std::vector<double> normalized_spectrum = normalizeVector(current_spectrum);
                std::vector<double> normalized_spectrum = current_spectrum;
                poi_value = calculateNRMSD(num_bins,ref_spectrum,normalized_spectrum);
            }
            else if (settings.parameter_of_interest == "chi_squared_g") {
                if (i_num == 0)
                    results_stream << settings.irradiation_conditions << ",";
                // std::vector<double> normalized_spectrum = normalizeVector(current_spectrum);
                std::vector<double> normalized_spectrum = current_spectrum;
                poi_value = calculateChiSquaredG(num_bins,ref_spectrum,normalized_spectrum);
            }
            else {
                throw std::logic_error("Unrecognized parameter of interest: " 
                    + settings.parameter_of_interest 
                    + ". Please refer to the README for allowed parameters"
                );
            }

            if (!settings.derivatives) {
                results_stream << poi_value;
                if (i_num == num_iteration_samples-1)
                    results_stream << "\n";

                else
                    results_stream << ",";
            }
            else {
                derivative_poi_values.push_back(poi_value);
            }
        }

        if (settings.derivatives) {
            std::vector<double> derivative_vector;
            calculateDerivatives(derivative_vector, num_iteration_samples, num_iterations_vector, derivative_poi_values);

            for (int i_num=0; i_num < num_iteration_samples; i_num++) {
                results_stream << derivative_vector[i_num];
                if (i_num == num_iteration_samples-1)
                    results_stream << "\n";
                else
                    results_stream << ",";
            }
        }

        // Save results for parameter of interest to CSV file
        std::ofstream output_file;
        output_file.open(settings.path_output_trend, std::ios_base::app);
        std::string results_string = results_stream.str();
        output_file << results_string;
        output_file.close();

        if (!settings.derivatives) {
            std::cout << "Saved 2D matrix of " << settings.parameter_of_interest << " values to " 
                << settings.path_output_trend << "\n";
        }
        else {
            std::cout << "Saved 2D matrix of derivatives of " << settings.parameter_of_interest 
                << " values to " << settings.path_output_trend << "\n";
        }

        // // Save results for parameter of interest to CSV file
        // std::ofstream output_file;
        // std::string settings.path_output_trend = settings.path_output_trend;
        // output_file.open(settings.path_output_trend, std::ios_base::app);
        // std::string results_string = results_stream.str();
        // output_file << results_string;
        // output_file.close();

        // std::cout << "Saved 2D matrix of " << settings.parameter_of_interest << " values to " << settings.path_output_trend << "\n";

    }

    //----------------------------------------------------------------------------------------------
    // Calculate some parameter of interest (POI) at specified numbers of iterations and beta value.
    // Iterate through a range of beta values, and a range of N values for each beta.
    // Output:
    // First row contains the iteration numbers
    // Other rows contain the beta value in the 1st column, followed by the POI value corresponding
    //  to the beta & N.
    // I.e. result is a 2D matrix of POI values (function of beta and N)
    // Visualize with plot_surface
    //----------------------------------------------------------------------------------------------
    if (settings.algorithm == "map") {
        // Create vector of beta values
        int num_orders_magnitude = log10(settings.beta_max/settings.beta_min);
        double current_beta = settings.beta_min;
        std::vector<double> beta_vector;
        for (int i=0; i<num_orders_magnitude; i++) {
            std::vector<double> temp_vector = linearSpacedDoubleVector(current_beta,current_beta*10,10);
            current_beta = current_beta*10;
            beta_vector.insert(beta_vector.end(), temp_vector.begin(), temp_vector.end());
        }

        // Create vector of number of iterations
        int num_increments = ((settings.iteration_max - settings.iteration_min) / settings.iteration_increment)+1;
        std::vector<int> num_iterations_vector = linearSpacedIntegerVector(
            settings.iteration_min,settings.iteration_max,num_increments
        );
        
        // Create needed variables
        int num_beta_samples = beta_vector.size();
        int num_iteration_samples = num_iterations_vector.size();

        std::vector<double> current_spectrum; // the reconstructed spectrum
        std::vector<double> energy_correction; // energy correction term used in MAP

        // Create stream to append results. First row is number of iteration increments
        std::ostringstream results_stream;
        // results_stream << std::scientific;
        results_stream << "0"; // empty first "cell"
        for (int i_num = 0; i_num < num_iteration_samples; i_num++) {
            results_stream << ",";
            results_stream << num_iterations_vector[i_num];
        }
        results_stream << "\n";

        // Loop through betas
        for (int i_beta=0; i_beta < num_beta_samples; i_beta++) {
            current_spectrum = initial_spectrum; // re-initialize spectrum for each beta
            results_stream << beta_vector[i_beta] << ",";

            // Loop through number of iterations
            for (int i_num=0; i_num < num_iteration_samples; i_num++) {
                int num_iterations;

                // only do # of iterations since previous run (i.e. don't do 3000, then 4000. Do
                // 3000 then iteration_size more, etc.)
                if (i_num == 0)
                    num_iterations = num_iterations_vector[i_num];
                else
                    num_iterations = num_iterations_vector[i_num]-num_iterations_vector[i_num-1];
                runMAP(energy_correction, beta_vector[i_beta], settings.prior, num_iterations, 
                    settings.error, num_measurements, num_bins, measurements, current_spectrum, 
                    nns_response, normalized_response, mlem_ratio
                );
            
                // Calculate one of the following parameters of interest & save to results stream
                double poi_value = 0;
                // Total fluence
                if (settings.parameter_of_interest == "total_fluence")
                    poi_value = calculateTotalFlux(num_bins,current_spectrum);
                // Total dose:
                else if (settings.parameter_of_interest == "total_dose")
                    poi_value = calculateDose(num_bins, current_spectrum, icrp_factors);
                // Total energy (penalty) term:
                else if (settings.parameter_of_interest == "total_energy_correction") {
                    poi_value = calculateTotalEnergyCorrection(energy_correction);
                }
                // Maximum "error" in MLEM ratio
                else if (settings.parameter_of_interest == "max_mlem_ratio") {
                    poi_value = calculateMaxRatio(num_measurements,mlem_ratio);
                }
                // Average "error" in MLEM ratio
                else if (settings.parameter_of_interest == "avg_mlem_ratio") {
                    poi_value = calculateAvgRatio(num_measurements,mlem_ratio);
                }
                else {
                    throw std::logic_error("Unrecognized parameter of interest: " 
                        + settings.parameter_of_interest 
                        + ". Please refer to the README for allowed parameters"
                    );
                }

                results_stream << poi_value;
                if (i_num == num_iteration_samples-1)
                    results_stream << "\n";

                else
                    results_stream << ",";
            }
        }

        // Save results for parameter of interest to CSV file
        std::ofstream output_file;
        output_file.open(settings.path_output_trend, std::ios_base::out);
        std::string results_string = results_stream.str();
        output_file << results_string;
        output_file.close();

        std::cout << "Saved 2D matrix of " << settings.parameter_of_interest << " values to " 
            << settings.path_output_trend << "\n";
    }

    return 0;
}



//**************************************************************************************************
// RETIRED/UNUSED ALGORITHM OPTIONS
//**************************************************************************************************

//----------------------------------------------------------------------------------------------
// For each spectral bin, determine the iteration index at which the MLEM correction factor is
// minimized.
// Outputs the iteration index for each spectral bin as a row
//----------------------------------------------------------------------------------------------
// if (settings.algorithm == "min_correction") {
//     // Create vector of number of iterations
//     int num_increments = ((settings.iteration_max - settings.iteration_min) / settings.iteration_increment)+1;
//     std::vector<int> num_iterations_vector = linearSpacedIntegerVector(
//        settings.iteration_min,settings.iteration_max,num_increments);
//     int num_iteration_samples = num_iterations_vector.size();

//     std::vector<double> current_spectrum = initial_spectrum; // the reconstructed spectrum

//     // Create stream to append results. First row is number of iteration increments
//     // determine if file exists
//     // std::string map_filename = "output/poi_output_mlem.csv";
//     std::ifstream rfile(settings.path_output_trend);
//     bool file_empty = is_empty(rfile);
//     // bool file_exists = rfile.good();
//     rfile.close();

//     // Add the energy bins to the file
//     std::ostringstream results_stream;
//     if (file_empty) {
//         results_stream << "Energy (MeV),";
//         for (int i_bin = 0; i_bin < num_bins; i_bin++) {
//             results_stream << energy_bins[i_bin];
//             if (i_bin != num_bins-1)
//                 results_stream << ",";
//         }
//         results_stream << "\n";
//     }

//     std::vector<double> min_deviations;
//     std::vector<int> min_indices;
//     double deviation = 0;

//     int total_num_iterations = 0;
//     for (int i_num=0; i_num < num_iteration_samples; i_num++) {
//         int num_iterations;

//         // only do # of iterations since previous run (i.e. don't do 3000, then 4000. Do
//         // 3000 then iteration_size more, etc.)
//         if (i_num == 0)
//             num_iterations = num_iterations_vector[i_num];
//         else
//             num_iterations = num_iterations_vector[i_num]-num_iterations_vector[i_num-1];
        // runMLEM(num_iterations, settings.error, num_measurements, num_bins, measurements, 
        //     current_spectrum, nns_response, normalized_response, mlem_ratio, mlem_correction, mlem_estimate
        // );

//         // on first iteration, save the deviations and iteration # for each bin
//         if (total_num_iterations == 0) {
//             for (int i_bin = 0; i_bin < num_bins; i_bin++) {
//                 min_deviations.push_back(abs(1.0 - mlem_correction[i_bin]));
//                 min_indices.push_back(num_iterations); 
//             }
//         }

//         total_num_iterations += num_iterations;
//         // on subsequent iterations, compare the current deviation with the minimum for each bin
//         // If new value is lower, overwrite existing minimum
//         // std::cout << "----------------------------------\n";
//         if (total_num_iterations > 0) {
//             for (int i_bin = 0; i_bin < num_bins; i_bin++) {
//                 deviation = abs(1.0 - mlem_correction[i_bin]);
//                 // std::cout << "i_bin=" << i_bin << ": dev=" << deviation << "\n";
//                 if (deviation < min_deviations[i_bin]) {
//                     min_deviations[i_bin] = deviation;
//                     min_indices[i_bin] = total_num_iterations; 
//                 }
//             }
//         }
//     }

//     results_stream << irradiation_conditions << ",";    
//     for (int i_bin = 0; i_bin < num_bins; i_bin++) {
//         results_stream << min_indices[i_bin];

//         if (i_bin != num_bins-1)
//             results_stream << ",";
//     }
//     results_stream << "\n";

//     // Save results for parameter of interest to CSV file
//     std::ofstream output_file;
//     std::string settings.path_output_trend = settings.path_output_trend;
//     output_file.open(settings.path_output_trend, std::ios_base::app);
//     std::string results_string = results_stream.str();
//     output_file << results_string;
//     output_file.close();

//     std::cout << "Saved indices of minimum deviations in correction factor to " << settings.path_output_trend << "\n";
// }



//----------------------------------------------------------------------------------------------
// Output the absolute difference between spectra obtained from subsequent iterations of
// unfolding (i.e. rows are num_bins elements long).
// Do this at specified iteration steps (i.e. one row for each sampled iteration step)
//----------------------------------------------------------------------------------------------
// if (settings.algorithm == "evolution") {
//     // Create vector of number of iterations
//     int num_increments = ((settings.iteration_max - settings.iteration_min) / settings.iteration_increment)+1;
//     std::vector<int> num_iterations_vector = linearSpacedIntegerVector(
//        settings.iteration_min,settings.iteration_max,num_increments);
//     int num_iteration_samples = num_iterations_vector.size();

//     std::vector<double> current_spectrum = initial_spectrum; // the reconstructed spectrum
//     std::vector<double> prev_spectrum;

//     // Add the energy bins to the file
//     std::ostringstream results_stream;
//     results_stream << "Energy (MeV),";
//     for (int i_bin = 0; i_bin < num_bins; i_bin++) {
//         results_stream << energy_bins[i_bin];
//         if (i_bin != num_bins-1)
//             results_stream << ",";
//     }
//     results_stream << "\n";

//     int total_num_iterations = 0;
//     for (int i_num=0; i_num < num_iteration_samples; i_num++) {
//         int num_iterations;

//         // only do # of iterations since previous run (i.e. don't do 3000, then 4000. Do
//         // 3000 then iteration_size more, etc.)
//         if (i_num == 0)
//             num_iterations = num_iterations_vector[i_num];
//         else
//             num_iterations = num_iterations_vector[i_num]-num_iterations_vector[i_num-1];
//         runMLEM_include_prev_spectrum(prev_spectrum,num_iterations, settings.error, 
        //     num_measurements, num_bins, measurements, current_spectrum, nns_response, normalized_response, 
        //     mlem_ratio, mlem_correction, mlem_estimate
        // );

//         total_num_iterations += num_iterations;

//         // Add the differences to the file
//         results_stream << "k = " << total_num_iterations << ",";
//         for (int i_bin = 0; i_bin < num_bins; i_bin++) {
//             results_stream << current_spectrum[i_bin] - prev_spectrum[i_bin];

//             if (i_bin != num_bins-1)
//                 results_stream << ",";
//         }
//         results_stream << "\n";
//     }        

//     // Save results for parameter of interest to CSV file
//     std::ofstream output_file;
//     std::string settings.path_output_trend = settings.path_output_trend;
//     output_file.open(settings.path_output_trend, std::ios_base::out);
//     std::string results_string = results_stream.str();
//     output_file << results_string;
//     output_file.close();

//     std::cout << "Saved spectral changes to " << settings.path_output_trend << "\n";
// }