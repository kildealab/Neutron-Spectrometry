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

    // Input filenames
    std::string config_file = "mlem.cfg";

    // Output filenames
    std::string dose_file = output_dir + "output_dose.csv";
    std::string o_spectrum_file = output_dir + "output_spectra.csv"; // result (unfolded) spectrum
    std::string report_file_pre = output_dir + "report_";
    std::string report_file_suf = ".txt";
    std::string figure_file_pre = output_dir + "figure_";
    std::string figure_file_suf = ".png";

    // Apply some settings read in from a config file
    int cutoff; // maximum # of MLEM itereations
    double norm; // vendor specfied normalization factor for the NNS used
    double error; // The target error on ratio between the experimental data points and the estimated data points from MLEM (e.g. 0.1 means the values must be within 10% of each other before the algorithm will terminate)
    double f_factor; // factor that converts measured charge (in fC) to counts per second [fA/cps]
    int num_poisson_samples;

    bool settings_success = setSettings(config_file, cutoff, norm, error, f_factor, num_poisson_samples);
    if (!settings_success) {
        throw std::logic_error("Unable to open configuration file: " + config_file);
    }
    double f_factor_report = f_factor; // original value read in
    f_factor = f_factor / 1e6; // Convert f_factor from fA/cps to nA/cps


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
    // Generate the inital spectrum matrix to input into MLEM algorithm:
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
    // Run the MLEM algorithm, iterating <cutoff> times.
    // Final result, i.e. MLEM-estimated spectrum, outputted in 'ini' matrix
    //----------------------------------------------------------------------------------------------
    std::vector<double> mlem_ratio; // vector that stores the ratio between measured data and MLEM estimated data
    int num_iterations = runMLEM(cutoff, error, num_measurements, num_bins, measurements, spectrum, nns_response, mlem_ratio);

    //----------------------------------------------------------------------------------------------
    // Display the result (output) matrix of the MLEM algorithm, which represents reconstructed 
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
    std::cout << "The ratios between measurements and MLEM-estimate measurements:" << '\n';

    for (int i_meas = 0; i_meas < num_measurements; i_meas++)
    {
        std::cout << mlem_ratio[i_meas] << '\n';
    }

    //----------------------------------------------------------------------------------------------
    // Display the number of iterations of MLEM that were actually executed (<= cutoff)
    //----------------------------------------------------------------------------------------------
    std::cout << '\n';
    std::cout << "The final number of MLEM iterations: " << num_iterations << std::endl;

    //----------------------------------------------------------------------------------------------
    // Repeat MLEM algorithm on Poisson-sampled dataset for num_poisson_samples iterations. Save the
    // sampled doses and spectra in arrays for subsequent RMSD uncertainty calcualtions.
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

        // Do MLEM on the initial spectrum & sampled measurement values
        runMLEM(cutoff, error, num_measurements, num_bins, sampled_measurements, sampled_spectrum, nns_response, sampled_mlem_ratio);
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
    // Calculate the ambient dose equivalent rate and its uncertainty
    //----------------------------------------------------------------------------------------------
    double ambient_dose_eq = calculateDose(num_bins, spectrum, icrp_factors);
    double ambient_dose_eq_uncertainty = calculateRMSD(num_poisson_samples, ambient_dose_eq, sampled_dose);

    //----------------------------------------------------------------------------------------------
    // Display the dose and its uncertainty
    //----------------------------------------------------------------------------------------------
    std::cout << '\n';
    std::cout << "The equivalent dose is: " << ambient_dose_eq << " mSv/h" << std::endl;
    std::cout << '\n';

    std::cout << '\n';
    std::cout << "The error on the equivalent dose is: " << ambient_dose_eq_uncertainty << " mSv/h" << std::endl;
    std::cout << '\n';

    //----------------------------------------------------------------------------------------------
    // Save results to file
    //----------------------------------------------------------------------------------------------
    saveDose(dose_file, irradiation_conditions, ambient_dose_eq, ambient_dose_eq_uncertainty);
    std::cout << "Saved calculated dose to " << dose_file << "\n";
    saveSpectrum(o_spectrum_file, irradiation_conditions, spectrum, spectrum_uncertainty, energy_bins);
    std::cout << "Saved unfolded spectrum to " << o_spectrum_file << "\n";

    std::string report_file = report_file_pre + irradiation_conditions + report_file_suf;
    prepareReport(report_file, irradiation_conditions, input_files, input_file_flags, cutoff, error, norm, f_factor_report, num_measurements, num_bins, num_poisson_samples, measurements_nc, dose_mu, doserate_mu, duration, energy_bins, initial_spectrum, nns_response, num_iterations, mlem_ratio, ambient_dose_eq, ambient_dose_eq_uncertainty, spectrum, spectrum_uncertainty, icrp_factors, GIT_COMMIT);
    std::cout << "Generated summary report: " << report_file << "\n\n";

    //----------------------------------------------------------------------------------------------
    // Plot the spectrum
    //----------------------------------------------------------------------------------------------
    plotSpectrum(figure_file_pre, figure_file_suf, irradiation_conditions, num_measurements, num_bins, energy_bins, spectrum, spectrum_uncertainty);

    return 0;
}