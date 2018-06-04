//**************************************************************************************************
// This program unfolds a set of measurements (in nC) obtained using a Nested Neutron Spectrometer
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

    // Names of input and output directories
    std::string input_dir = "input/";
    std::string output_dir = "output/";

    // NOTE: Indices are linked between the following arrays and vectors (i.e. input_files[0]
    // corresponds to input_file_flags[0] and input_file_defaults[0])
    // Array that stores the allowed options that specify input files
    // Add new options at end of array
    const int num_ifiles = 5;
    std::string input_file_flags_arr[num_ifiles] = {
        "--measurements",
        "--input-spectrum",
        "--energy-bins",
        "--nns-response",
        "--icrp-factors"
    };
    // Array that stores default filename for each input file
    std::string input_file_defaults_arr[num_ifiles] = {
        "measurements.txt",
        "spectrum_step.csv",
        "energy_bins.csv",
        "nns_response.csv",
        "icrp_conversions.csv"
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
        setfile(arg_vector, input_dir, input_file_flags[i], input_file_defaults[i], input_files[i]);
    }

    // Notify user if unknown parameters were received
    checkUnknownParameters(arg_vector, input_file_flags);

    // Apply some settings read in from a config file
    UnfoldingSettings settings;
    setSettings(input_dir + "auto.cfg", settings);

    settings.set_f_factor(settings.f_factor / 1e6); // Convert f_factor from fA/cps to nA/cps

    // Standard figure file suffix (file extension)
    std::string figure_file_suf = ".png";

    // Read measured data (in nC) from input file
    std::string irradiation_conditions;
    double dose_mu; // dose delivered (MU) for individual measurement
    double doserate_mu; // dose rate (MU/min) used for individual measurement
    int duration; // Duration (s) of individual measurement acquisition
    std::vector<double> measurements_nc = getMeasurements(input_files[0], irradiation_conditions, dose_mu, doserate_mu, duration);
    int num_measurements = measurements_nc.size();

    // Convert measured charge in nC to counts per second
    // Re-order measurments from (7 moderators to 0) to (0 moderators to 7)
    std::vector<double> measurements;
    for (int index=0; index < num_measurements; index++) {
        double measurement_cps = measurements_nc[num_measurements-index-1]*settings.norm/settings.f_factor/duration;
        measurements.push_back(measurement_cps);
    }

    //----------------------------------------------------------------------------------------------
    // Generate the energy bins matrix:
    //  - size = # of energy bins
    // Input the energies from energy bins file: input_files[2]
    //  - values in units of [MeV]
    //----------------------------------------------------------------------------------------------
    std::vector<double> energy_bins;
    readInputFile1D(input_files[2],energy_bins);

    int num_bins = energy_bins.size();

    //----------------------------------------------------------------------------------------------
    // Generate the detector response matrix (representing the dector response function):
    //  - outer size = # of measurements
    //  - inner size = # of energy bins
    // Detector response value for each # of moderators for each energy (currently 52 energies)
    // Input the response from input_files[3]
    //  - values in units of [cm^2]
    //
    // The response function accounts for variable number of (n,p) reactions in He-3 for each
    // moderators, as a function of energy. Calculated by vendor using MC
    //----------------------------------------------------------------------------------------------
    std::vector<std::vector<double>> nns_response;
    readInputFile2D(input_files[3],nns_response);
    checkDimensions(num_measurements, "number of measurements", nns_response.size(), "NNS response");
    checkDimensions(num_bins, "number of energy bins", nns_response[0].size(), "NNS response");

    //----------------------------------------------------------------------------------------------
    // Generate the inital spectrum matrix to input into the unfolding algorithm:
    //  - size = # of energy bins
    // Input from the input spectrum file (input_files[1])
    //  - values are neutron fluence rates [neutrons cm^-2 s^-1])
    //  - Currently (2017-08-16) input a step function (high at thermals & lower), because a flat 
    //  spectrum underestimates (does not yield any) thermal neutrons
    //----------------------------------------------------------------------------------------------
    std::vector<double> initial_spectrum;
    readInputFile1D(input_files[1],initial_spectrum);
    checkDimensions(num_bins, "number of energy bins", initial_spectrum.size(), "Input spectrum");

    std::vector<double> spectrum = initial_spectrum; // save the initial spectrum for report output

    //----------------------------------------------------------------------------------------------
    // Generate the ICRP conversion matrix (factors to convert fluence to ambient dose equivalent):
    //  - size = # of energy bins
    // Input from input_files[4]
    //  - values are in units of [pSv cm^2]
    //  - H values were obtained by linearly interopolating tabulated data to match energy bins used
    // Page 200 of document (ICRP 74 - ATables.pdf)
    //----------------------------------------------------------------------------------------------
    std::vector<double> icrp_factors;
    readInputFile1D(input_files[4],icrp_factors);
    checkDimensions(num_bins, "number of energy bins", icrp_factors.size(), "Number of ICRP factors");

    //----------------------------------------------------------------------------------------------
    // Run the automatic unfolding algorithm.
    // Note: the normalized system matrix is calculated first. It is required in unfolding, and is
    // a constant value.
    //----------------------------------------------------------------------------------------------
    std::vector<double> normalized_response = normalizeResponse(num_bins, num_measurements, nns_response);

    std::vector<double> mlem_ratio; // vector that stores the ratio between measured data and MLEM estimated data

    // If doing MLEM unfolding:
    if (settings.algorithm == "mlem") {
        // Create vector of number of iterations
        int num_increments = ((settings.max_num_iterations - settings.min_num_iterations) / settings.iteration_increment)+1;
        std::vector<int> num_iterations_vector = linearSpacedIntegerVector(settings.min_num_iterations,settings.max_num_iterations,num_increments);
        int num_iteration_samples = num_iterations_vector.size();

        std::vector<double> current_spectrum = initial_spectrum; // the reconstructed spectrum

        // Create stream to append results. First row is number of iteration increments
        // determine if file exists
        std::string map_filename = "output/poi_output_mlem.csv";
        std::ifstream rfile(map_filename);
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

        // Loop through number of iterations
        for (int i_num=0; i_num < num_iteration_samples; i_num++) {
            int num_iterations;

            // only do # of iterations since previous run (i.e. don't do 3000, then 4000. Do
            // 3000 then iteration_size more, etc.)
            if (i_num == 0)
                num_iterations = num_iterations_vector[i_num];
            else
                num_iterations = num_iterations_vector[i_num]-num_iterations_vector[i_num-1];
            runMLEM(num_iterations, settings.error, num_measurements, num_bins, measurements, current_spectrum, nns_response, normalized_response, mlem_ratio);
        
            // Calculate one of the following parameters of interest & save to results stream
            double poi_value = 0;
            // Total fluence
            if (settings.parameter_of_interest == "total_fluence") {
                if (i_num == 0)
                    results_stream << "Total fluence,";
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
                    results_stream << "Maximum MLEM ratio,";
                poi_value = calculateMaxRatio(num_measurements,mlem_ratio);
            }
            // Average "error" in MLEM ratio
            else if (settings.parameter_of_interest == "avg_mlem_ratio") {
                if (i_num == 0)
                    results_stream << "Average MLEM ratio,";
                poi_value = calculateAvgRatio(num_measurements,mlem_ratio);
            }
            else {
                throw std::logic_error("Unrecognized parameter of interest: " + settings.parameter_of_interest + ". Please refer to the README for allowed parameters");
            }

            results_stream << poi_value;
            if (i_num == num_iteration_samples-1)
                results_stream << "\n";

            else
                results_stream << ",";
        }

        // Save results for parameter of interest to CSV file
        std::ofstream map_file;
        map_file.open(map_filename, std::ios_base::app);
        std::string results_string = results_stream.str();
        map_file << results_string;
        map_file.close();

        std::cout << "Saved 2D matrix of " << settings.parameter_of_interest << " values to " << map_filename << "\n";
    }
    // If doing MAP unfolding
    if (settings.algorithm == "map") {
        // Create vector of beta values
        int num_orders_magnitude = log10(settings.max_beta/settings.min_beta);
        double current_beta = settings.min_beta;
        std::vector<double> beta_vector;
        for (int i=0; i<num_orders_magnitude; i++) {
            std::vector<double> temp_vector = linearSpacedDoubleVector(current_beta,current_beta*10,10);
            current_beta = current_beta*10;
            beta_vector.insert(beta_vector.end(), temp_vector.begin(), temp_vector.end());
        }

        // Create vector of number of iterations
        int num_increments = ((settings.max_num_iterations - settings.min_num_iterations) / settings.iteration_increment)+1;
        std::vector<int> num_iterations_vector = linearSpacedIntegerVector(settings.min_num_iterations,settings.max_num_iterations,num_increments);
        
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
                runMAP(energy_correction, beta_vector[i_beta], settings.prior, num_iterations, settings.error, num_measurements, num_bins, measurements, current_spectrum, nns_response, normalized_response, mlem_ratio);
            
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
                    throw std::logic_error("Unrecognized parameter of interest: " + settings.parameter_of_interest + ". Please refer to the README for allowed parameters");
                }

                results_stream << poi_value;
                if (i_num == num_iteration_samples-1)
                    results_stream << "\n";

                else
                    results_stream << ",";
            }
        }

        // Save results for parameter of interest to CSV file
        std::string map_filename = "output/poi_output_map.csv";
        std::ofstream map_file;
        map_file.open(map_filename, std::ios_base::out);
        std::string results_string = results_stream.str();
        map_file << results_string;
        map_file.close();

        std::cout << "Saved 2D matrix of " << settings.parameter_of_interest << " values to " << map_filename << "\n";
    }

    return 0;
}