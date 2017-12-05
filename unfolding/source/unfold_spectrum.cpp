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
    // Remove file extension and director
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
    int cutoff; // maximum # of iterations in the unfolding algorithm
    double norm; // vendor specfied normalization factor for the NNS used
    double error; // The target error on ratio between the experimental data points and the estimated data points from unfolding (e.g. 0.1 means the values must be within 10% of each other before the algorithm will terminate)
    double f_factor; // factor that converts measured charge (in fC) to counts per second [fA/cps]
    double beta; // factor that affects damping of high noise spectral component in MAP method
    int num_poisson_samples;

    bool settings_success = setSettings(input_files[5], algorithm_name, cutoff, norm, error, f_factor, beta, num_poisson_samples);
    if (!settings_success) {
        throw std::logic_error("Unable to open configuration file: " + input_files[5]);
    }
    double f_factor_report = f_factor; // original value read in
    f_factor = f_factor / 1e6; // Convert f_factor from fA/cps to nA/cps

    // Standard figure file suffix (file extension)
    std::string figure_file_suf = ".png";
    // Include algorithm name and cutoff # of iterations in the figure file suffix
    // std::ostringstream temp_stream;
    // temp_stream << algorithm_name << cutoff << ".png";
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
        double measurement_cps = measurements_nc[num_measurements-index-1]*norm/f_factor/duration;
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
        num_iterations = runMLEM(cutoff, error, num_measurements, num_bins, measurements, spectrum, nns_response, normalized_response, mlem_ratio);
    }
    else if (algorithm_name == "map") {
        std::vector<double> energy_correction;
        num_iterations = runMAP(energy_correction, beta, cutoff, error, num_measurements, num_bins, measurements, spectrum, nns_response, normalized_response, mlem_ratio);
    }
    // else if (algorithm_name == "auto") {
    //     std::vector<double> energy_correction;

    //     static const double beta_array[] = {1e-11,5e-11,1e-10,5e-10,1e-9,5e-9,1e-8,5e-8,1e-7};
    //     static const int n_array[] = {1000,2000,3000,4000,5000,6000,7000,8000};

    //     int betas_length = sizeof(beta_array)/sizeof(*beta_array);
    //     int n_length = sizeof(n_array)/sizeof(*n_array);

    //     std::vector<double> test_spectrum;
    //     bool found_beta = false;

    //     // Find beta value
    //     for (int i_beta=0; i_beta < betas_length; i_beta++) {
    //         test_spectrum = initial_spectrum; // reset spectrum on each beta
    //         for (int i_n=0; i_n < n_length; i_n++) {
    //             // only do # of iterations since previous run (i.e. don't do 3000, then 4000. Do
    //             // 3000 then 1000 more, etc.)
    //             int n_iterations;
    //             if (i_n == 0)
    //                 n_iterations = n_array[i_n];
    //             else
    //                 n_iterations = n_array[i_n]-n_array[i_n-1];
    //             runMAP(energy_correction, beta_array[i_beta], n_array[i_n], error, num_measurements, num_bins, measurements, test_spectrum, nns_response, normalized_response, mlem_ratio);
                
    //             // get max energy correction term
    //             double beta_threshold = 1e-3;
    //             double max_correction = 0.0;
    //             for (int i_bin = 0; i_bin < num_bins; i_bin++)
    //             {
    //                 if (energy_correction[i_bin]>max_correction) {
    //                     max_correction = energy_correction[i_bin];
    //                 }
    //                 if (energy_correction[i_bin]>beta_threshold) {
    //                     beta = beta_array[i_beta];
    //                     found_beta = true;
    //                     break;
    //                 }
    //             }
    //             std::cout << "cur_beta: " << beta_array[i_beta] << " | cur_n: " << n_array[i_n] << " | max correction: " << max_correction << "\n";
    //             if (found_beta)
    //                 break;
    //         }
    //         if (found_beta)
    //             break;
    //     }

    //     std::cout << "Found beta? " << found_beta << "\n";
    //     std::cout << "Beta " << beta << "\n";
    //     return 0;
    // }
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

    for (int i_poiss = 0; i_poiss < num_poisson_samples; i_poiss++) {
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
            runMLEM(cutoff, error, num_measurements, num_bins, sampled_measurements, sampled_spectrum, nns_response, normalized_response, sampled_mlem_ratio);
        }
        else if (algorithm_name == "map") {
            std::vector<double> sampled_energy_correction;
            runMAP(sampled_energy_correction, beta, cutoff, error, num_measurements, num_bins, sampled_measurements, sampled_spectrum, nns_response, normalized_response, sampled_mlem_ratio);
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
    calculateRMSD_vector(num_poisson_samples, spectrum, sampled_spectra, spectrum_uncertainty);

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
    double ambient_dose_eq_uncertainty = calculateRMSD(num_poisson_samples, ambient_dose_eq, sampled_dose);

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
    prepareReport(report_file, irradiation_conditions, input_files, input_file_flags, algorithm_name, cutoff, error, norm, f_factor_report, beta, num_measurements, num_bins, num_poisson_samples, measurements_nc, dose_mu, doserate_mu, duration, energy_bins, initial_spectrum, nns_response, num_iterations, mlem_ratio, ambient_dose_eq, ambient_dose_eq_uncertainty, spectrum, spectrum_uncertainty, icrp_factors, GIT_COMMIT);
    std::cout << "Generated summary report: " << report_file << "\n\n";

    //----------------------------------------------------------------------------------------------
    // Plot the spectrum
    //----------------------------------------------------------------------------------------------
    plotSpectrum(figure_file_pre, figure_file_suf, irradiation_conditions, num_measurements, num_bins, energy_bins, spectrum, spectrum_uncertainty);

    return 0;
}