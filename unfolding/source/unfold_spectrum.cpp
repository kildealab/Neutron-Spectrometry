//**************************************************************************************************
// This program unfolds a set of measurements (in nC) obtained using a Nested Neutron Spectrometer
// into a spectrum of neutron flux. Also,the program calculates the ambient dose equivalent rate
// associated with the unfolded spectrum. Uncertainty estimates for the spectrum and dose are
// provided. Input parameters (measurements, unfolding algorithm configuration, etc.) are read in
// from files located in the "input/" directory. Results are output into the "output/" directory,
// and include:
//  - the dose and its uncertainty
//  - the spectrum and its uncertainty (numeric and graphical forms)
//  - a report that details the execution of this program for archival and reproducibility
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
    const int num_ifiles = 6;
    std::string input_file_flags_arr[num_ifiles] = {
        "--measurements",
        "--input-spectrum",
        "--energy-bins",
        "--nns-response",
        "--icrp-factors",
        "--configuration"
    };
    // Array that stores default filename for each input file
    std::string input_file_defaults_arr[num_ifiles] = {
        "measurements.txt",
        "spectrum_step.csv",
        "energy_bins.csv",
        "nns_response.csv",
        "icrp_conversions.csv",
        "mlem.cfg"
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

    // Get the name of the unfolding algorithm
    // Remove file extension and directory
    std::string algorithm_name = input_files[5].substr(0, input_files[5].find(".", 0));
    algorithm_name = algorithm_name.substr(algorithm_name.find("/", 0)+1);
    transform(algorithm_name.begin(), algorithm_name.end(), algorithm_name.begin(),(int (*)(int))tolower);

    // Output filenames
    std::string dose_file = output_dir + "output_dose.csv";
    std::string o_spectrum_file = output_dir + "output_spectra.csv"; // result (unfolded) spectrum
    std::string report_file_pre = output_dir + "report_";
    std::string report_file_suf = ".txt";
    std::string figure_file_pre = output_dir + "figure_";

    // Apply some settings read in from a config file
    UnfoldingSettings settings;
    setSettings(input_files[5], algorithm_name, settings);

    double f_factor_report = settings.f_factor; // original value read in
    settings.set_f_factor(settings.f_factor / 1e6); // Convert f_factor from fA/cps to nA/cps

    // Standard figure file suffix (file extension)
    std::string figure_file_suf = ".png";
    // Include algorithm name and cutoff # of iterations in the figure file suffix
    // std::ostringstream temp_stream;
    // temp_stream << algorithm_name << settings.cutoff << ".png";
    // std::string figure_file_suf = temp_stream.str();

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
    // Print out the processed measured data matrix
    //----------------------------------------------------------------------------------------------ls
    std::cout << '\n';
    std::cout << "The measurements in CPS are:" << '\n'; // newline

    //Loop over the data matrix, display each value
    for (int i = 0; i < 8; ++i)
    {
        std::cout << measurements[i] << '\n';
    }
    std::cout << '\n';

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
    // Run the unfolding algorithm, iterating <cutoff> times.
    // Final result, i.e. unfolded spectrum, outputted in 'ini' matrix
    // Note: the normalized system matrix is calculated first. It is required in unfolding, and is
    // a constant value.
    //----------------------------------------------------------------------------------------------
    std::vector<double> normalized_response = normalizeResponse(num_bins, num_measurements, nns_response);

    std::vector<double> mlem_ratio; // vector that stores the ratio between measured data and MLEM estimated data
    int num_iterations;

    if (algorithm_name == "mlem") {
        num_iterations = runMLEM(settings.cutoff, settings.error, num_measurements, num_bins, measurements, spectrum, nns_response, normalized_response, mlem_ratio);
    }
    else if (algorithm_name == "map") {
        std::vector<double> energy_correction;
        num_iterations = runMAP(energy_correction, settings.beta, settings.prior, settings.cutoff, settings.error, num_measurements, num_bins, measurements, spectrum, nns_response, normalized_response, mlem_ratio);
        // return 1;
    }
    else if (algorithm_name == "auto") {
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
                    for (auto& n : energy_correction)
                        poi_value += n;
                }
                // Maximum "error" in MLEM ratio
                else if (settings.parameter_of_interest == "max_mlem_ratio") {
                    double max_mlem_ratio = 0.0;
                    for (int i=0; i < num_measurements; i++) {
                        if (mlem_ratio[i] > max_mlem_ratio) {
                            max_mlem_ratio = mlem_ratio[i];
                        }
                    }
                    poi_value = max_mlem_ratio;
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
        std::string map_filename = "output/map_file.csv";
        std::ofstream map_file;
        map_file.open(map_filename, std::ios_base::out);
        std::string results_string = results_stream.str();
        map_file << results_string;
        map_file.close();

        std::cout << "Saved 2D matrix of " << settings.parameter_of_interest << " values to " << map_filename << "\n";

        return 1;
    }
    else {
        //throw error
        std::cout << "No unfolding algorithm found for: " + algorithm_name + '\n';
    }

    //----------------------------------------------------------------------------------------------
    // Display the result (output) matrix of the unfolding algorithm, which represents reconstructed 
    // spectral data.
    //----------------------------------------------------------------------------------------------
    std::cout << "The unfolded spectrum:" << '\n';

    for (int i_bin = 0; i_bin < num_bins; i_bin++)
    {
        std::cout << spectrum[i_bin] << '\n';
    }

    //----------------------------------------------------------------------------------------------
    // Display the ratio (error) matrix of the MLEM algorithm, which represents the deviation of the
    // MLEM-generated measured charge values (which correspond to the above MLEM-generated measured 
    // spectrum) from the actual measured charge values. 
    //----------------------------------------------------------------------------------------------
    std::cout << '\n';
    std::cout << "The ratios between measurements and unfolded measurements:" << '\n';

    for (int i_meas = 0; i_meas < num_measurements; i_meas++)
    {
        std::cout << mlem_ratio[i_meas] << '\n';
    }

    //----------------------------------------------------------------------------------------------
    // Display the number of unfolding iterations that were actually executed (<= cutoff)
    //----------------------------------------------------------------------------------------------
    std::cout << '\n';
    std::cout << "The final number of unfolding iterations: " << num_iterations << std::endl;

    //----------------------------------------------------------------------------------------------
    // Repeat unfolding algorithm on Poisson-sampled dataset for num_poisson_samples iterations.
    // Save the sampled doses and spectra in arrays for subsequent RMSD uncertainty calcualtions.
    //----------------------------------------------------------------------------------------------
    std::vector<std::vector<double>> sampled_spectra; // dimensions: num_poisson_samples x num_bins
    std::vector<double> sampled_dose; // dimension: num_poisson_samples

    for (int i_poiss = 0; i_poiss < settings.num_poisson_samples; i_poiss++) {
        std::vector<double> sampled_measurements; // dimension: num_measurements
        std::vector<double> sampled_mlem_ratio; // dimension: num_measurements
        std::vector<double> sampled_spectrum = initial_spectrum; // dimension: num_bins

        // Create Poisson sampled measurement values (CPS)
        for (int i_meas = 0; i_meas < num_measurements; i_meas++) {
            double sampled_value = poisson(measurements[i_meas]);
            sampled_measurements.push_back(sampled_value);
        }

        // Do unfolding on the initial spectrum & sampled measurement values
        if (algorithm_name == "mlem") {
            runMLEM(settings.cutoff, settings.error, num_measurements, num_bins, sampled_measurements, sampled_spectrum, nns_response, normalized_response, sampled_mlem_ratio);
        }
        else if (algorithm_name == "map") {
            std::vector<double> sampled_energy_correction;
            runMAP(sampled_energy_correction, settings.beta, settings.prior, settings.cutoff, settings.error, num_measurements, num_bins, sampled_measurements, sampled_spectrum, nns_response, normalized_response, sampled_mlem_ratio);
        }
        else {
            //throw error
            std::cout << "No unfolding algorithm found for: " + algorithm_name + '\n';
        }        

        sampled_spectra.push_back(sampled_spectrum); // add to growing array of sampled spectra

        // Calculate the ambient dose equivalent associated with the sampled spectrum
        double sdose = calculateDose(num_bins, sampled_spectrum, icrp_factors);

        sampled_dose.push_back(sdose);
    }

    //----------------------------------------------------------------------------------------------
    // Calculate RMSD between sampled spectra and the unfolded spectrum
    //----------------------------------------------------------------------------------------------
    std::vector<double> spectrum_uncertainty;
    calculateRMSD_vector(settings.num_poisson_samples, spectrum, sampled_spectra, spectrum_uncertainty);

    //----------------------------------------------------------------------------------------------
    // Print the uncertainty spectrum
    //----------------------------------------------------------------------------------------------
    // std::cout << '\n'; // newline
    // std::cout << "The uncertainy in the measured spectrum:" << '\n'; // newline

    // for (int i = 0; i < 52; i++)
    // {
    //     std::cout << spectrum_uncertainty[i] << '\n';
    // }

    //----------------------------------------------------------------------------------------------
    // Calculate quantities of interest (e.g. dose & its uncertainty)
    //----------------------------------------------------------------------------------------------
    double ambient_dose_eq = calculateDose(num_bins, spectrum, icrp_factors);
    double ambient_dose_eq_uncertainty = calculateRMSD(settings.num_poisson_samples, ambient_dose_eq, sampled_dose);

    double total_charge = calculateTotalCharge(num_measurements,measurements_nc);

    double total_flux = calculateTotalFlux(num_bins,spectrum);
    double total_flux_uncertainty = calculateSumUncertainty(num_bins,spectrum_uncertainty);

    double avg_energy = calculateAverageEnergy(num_bins,spectrum,energy_bins);
    double avg_energy_uncertainty = calculateEnergyUncertainty(num_bins,energy_bins,spectrum,spectrum_uncertainty,total_flux,total_flux_uncertainty);

    // double source_strength = calculateSourceStrength(num_bins,spectrum,duration,dose_mu);

    //----------------------------------------------------------------------------------------------
    // Display calculated quantities
    //----------------------------------------------------------------------------------------------
    std::cout << '\n';
    std::cout << "The equivalent dose is: " << ambient_dose_eq << " mSv/h" << std::endl;
    std::cout << "The uncertainty on the equivalent dose is: " << ambient_dose_eq_uncertainty << " mSv/h" << std::endl;
    std::cout << '\n';

    std::cout << "The total measured charge is: " << total_charge << " nC" << std::endl;
    std::cout << '\n';

    std::cout << "The total neutron flux is: " << total_flux << " n cm^-2 s^-1" << std::endl;
    std::cout << "The uncertainty on the total flux is: " << total_flux_uncertainty << " n cm^-2 s^-1" << std::endl;
    std::cout << '\n';

    std::cout << "The average neutron energy is: " << avg_energy << " MeV" << std::endl;
    std::cout << "The uncertainty on the average energy is: " << avg_energy_uncertainty << " MeV" << std::endl;
    std::cout << '\n';

    // std::cout << "The neutron source strength is: " << source_strength << " n Gy^-1" << std::endl;
    // std::cout << '\n';


    //----------------------------------------------------------------------------------------------
    // Save results to file
    //----------------------------------------------------------------------------------------------
    saveDose(dose_file, irradiation_conditions, ambient_dose_eq, ambient_dose_eq_uncertainty);
    std::cout << "Saved calculated dose to " << dose_file << "\n";
    saveSpectrumAsRow(o_spectrum_file, num_bins, irradiation_conditions, spectrum, spectrum_uncertainty, energy_bins);
    std::cout << "Saved unfolded spectrum to " << o_spectrum_file << "\n";

    std::string report_file = report_file_pre + irradiation_conditions + report_file_suf;
    prepareReport(report_file, irradiation_conditions, input_files, input_file_flags, algorithm_name, settings.cutoff, settings.error, settings.norm, f_factor_report, settings.beta, num_measurements, num_bins, settings.num_poisson_samples, measurements_nc, dose_mu, doserate_mu, duration, energy_bins, initial_spectrum, nns_response, num_iterations, mlem_ratio, ambient_dose_eq, ambient_dose_eq_uncertainty, total_charge, total_flux, total_flux_uncertainty, avg_energy, avg_energy_uncertainty, spectrum, spectrum_uncertainty, icrp_factors, GIT_COMMIT);
    std::cout << "Generated summary report: " << report_file << "\n\n";

    //----------------------------------------------------------------------------------------------
    // Plot the spectrum
    //----------------------------------------------------------------------------------------------
    plotSpectrum(figure_file_pre, figure_file_suf, irradiation_conditions, num_measurements, num_bins, energy_bins, spectrum, spectrum_uncertainty);

    return 0;
}