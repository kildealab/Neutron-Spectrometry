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
    const int num_ifiles = 1;
    std::string input_file_flags_arr[num_ifiles] = {
        "--configuration"
    };
    // Array that stores default filename for each input file
    std::string input_file_defaults_arr[num_ifiles] = {
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

    // Output filenames    
    std::string figure_file_pre = output_dir + "figure_";

    // Apply some settings read in from a config file
    UnfoldingSettings settings;
    setSettings(input_files[0], settings);

    double f_factor_report = settings.f_factor; // original value read in
    settings.set_f_factor(settings.f_factor / 1e6); // Convert f_factor from fA/cps to nA/cps

    // Standard figure file suffix (file extension)
    std::string figure_file_suf = ".png";
    // Include algorithm name and cutoff # of iterations in the figure file suffix
    // std::ostringstream temp_stream;
    // temp_stream << algorithm_name << settings.cutoff << ".png";
    // std::string figure_file_suf = temp_stream.str();

    // Read measured data from input file
    std::string irradiation_conditions;
    double dose_mu; // dose delivered (MU) for individual measurement
    double doserate_mu; // dose rate (MU/min) used for individual measurement
    int duration; // Duration (s) of individual measurement acquisition

    std::vector<double> measurements_nc;
    std::vector<double> measurements;
    int num_measurements = 0;

    // Handle inputs with different units
    if (settings.meas_units == "cps") {
        measurements = getMeasurementsCPS(settings.measurements_path, irradiation_conditions);
        num_measurements = measurements.size();

        // for (int index=0; index < num_measurements; index++) {
        //     measurements[index] = measurements[index]*settings.norm;
        // }
        std::reverse(measurements.begin(),measurements.end());
    }
    else {
        measurements_nc = getMeasurements(settings.measurements_path, irradiation_conditions, dose_mu, doserate_mu, duration);
        num_measurements = measurements_nc.size();

        for (int index=0; index < num_measurements; index++) {
            double measurement_cps = measurements_nc[num_measurements-index-1]*settings.norm/settings.f_factor/duration;
            measurements.push_back(measurement_cps);
        }
    }

    //----------------------------------------------------------------------------------------------
    // Print out the processed measured data matrix
    //----------------------------------------------------------------------------------------------
    std::cout << '\n';
    std::cout << "The measurements in CPS are:" << '\n'; // newline

    //Loop over the data matrix, display each value
    for (int i = 0; i < num_measurements; ++i)
    {
        std::cout << measurements[i] << '\n';
    }
    std::cout << '\n';

    //----------------------------------------------------------------------------------------------
    // Generate the energy bins matrix:
    //  - size = # of energy bins
    // Input the energies from energy bins file
    //  - values in units of [MeV]
    //----------------------------------------------------------------------------------------------
    std::vector<double> energy_bins;
    readInputFile1D(settings.energy_bins_path,energy_bins);

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
    readInputFile2D(settings.system_response_path,nns_response);
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
    readInputFile1D(settings.input_spectrum_path,initial_spectrum);
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
    readInputFile1D(settings.icrp_factors_path,icrp_factors);
    checkDimensions(num_bins, "number of energy bins", icrp_factors.size(), "Number of ICRP factors");

    //----------------------------------------------------------------------------------------------
    // Run the unfolding algorithm, iterating <cutoff> times.
    // Final result, i.e. unfolded spectrum, outputted in 'ini' matrix
    // Note: the normalized system matrix is calculated first. It is required in unfolding, and is
    // a constant value.
    //----------------------------------------------------------------------------------------------
    std::vector<double> normalized_response = normalizeResponse(num_bins, num_measurements, nns_response);

    std::vector<double> mlem_ratio; // vector that stores the ratio between measured data and MLEM estimated data
    std::vector<double> mlem_correction; // vector that stores the correction factors applied in each spectral bin
    std::vector<double> mlem_estimate; // vector that stores the MLEM estimated data
    int num_iterations;

    // MLEM-STOP parameters
    double j_factor = 0;
    double j_threshold = 0;

    if (settings.algorithm == "mlem") {
        num_iterations = runMLEM(settings.cutoff, settings.error, num_measurements, num_bins,
            measurements, spectrum, nns_response, normalized_response, mlem_ratio, mlem_correction, 
            mlem_estimate
        );
    }
    else if (settings.algorithm == "mlemstop") {
        j_threshold = determineJThreshold(num_measurements,measurements,settings.cps_crossover);

        num_iterations = runMLEMSTOP(settings.cutoff, num_measurements, num_bins, measurements,
            spectrum, nns_response, normalized_response, mlem_ratio, mlem_correction, mlem_estimate,
            j_threshold, j_factor
        );
    }
    else if (settings.algorithm == "map") {
        std::vector<double> energy_correction;
        num_iterations = runMAP(energy_correction, settings.beta, settings.prior, settings.cutoff, 
            settings.error, num_measurements, num_bins, measurements, spectrum, nns_response, 
            normalized_response, mlem_ratio
        );
    }
    else {
        //throw error
        std::cout << "No unfolding algorithm found for: " + settings.algorithm + '\n';
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
    std::cout << "The reconstructed measured values (CPS):" << '\n';

    for (int i_meas = 0; i_meas < num_measurements; i_meas++)
    {
        std::cout << mlem_estimate[i_meas] << '\n';
    }

    std::cout << '\n';
    std::cout << "The ratios between measurements and reconstructed measurements:" << '\n';

    for (int i_meas = 0; i_meas < num_measurements; i_meas++)
    {
        std::cout << mlem_ratio[i_meas] << '\n';
        // std::cout << abs(1.0 - mlem_ratio[i_meas]) << '\n';
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
        std::vector<double> sampled_mlem_correction; // dimension: num_measurements
        std::vector<double> sampled_mlem_estimate; // dimension: num_measurements
        std::vector<double> sampled_spectrum = initial_spectrum; // dimension: num_bins

        // Create Poisson sampled measurement values (CPS)
        for (int i_meas = 0; i_meas < num_measurements; i_meas++) {
            double sampled_value = poisson(measurements[i_meas]);
            sampled_measurements.push_back(sampled_value);
        }

        // Do unfolding on the initial spectrum & sampled measurement values
        if (settings.algorithm == "mlem" || settings.algorithm == "mlemstop") {
            runMLEM(num_iterations, settings.error, num_measurements, num_bins, sampled_measurements, 
                sampled_spectrum, nns_response, normalized_response, sampled_mlem_ratio, 
                sampled_mlem_correction, sampled_mlem_estimate
            );
        }
        else if (settings.algorithm == "map") {
            std::vector<double> sampled_energy_correction;
            runMAP(sampled_energy_correction, settings.beta, settings.prior, settings.cutoff, 
                settings.error, num_measurements, num_bins, sampled_measurements, sampled_spectrum, 
                nns_response, normalized_response, sampled_mlem_ratio
            );
        }
        else {
            //throw error
            std::cout << "No unfolding algorithm found for: " + settings.algorithm + '\n';
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

    // double total_charge = calculateTotalCharge(num_measurements,measurements_nc);

    double total_flux = calculateTotalFlux(num_bins,spectrum);
    double total_flux_uncertainty = calculateSumUncertainty(num_bins,spectrum_uncertainty);

    double avg_energy = calculateAverageEnergy(num_bins,spectrum,energy_bins);
    double avg_energy_uncertainty = calculateEnergyUncertainty(num_bins,energy_bins,spectrum,
        spectrum_uncertainty,total_flux,total_flux_uncertainty
    );

    // double source_strength = calculateSourceStrength(num_bins,spectrum,duration,dose_mu);

    //----------------------------------------------------------------------------------------------
    // Display calculated quantities
    //----------------------------------------------------------------------------------------------
    std::cout << '\n';
    std::cout << "The equivalent dose is: " << ambient_dose_eq << " mSv/h" << std::endl;
    std::cout << "The uncertainty on the equivalent dose is: " << ambient_dose_eq_uncertainty << " mSv/h" << std::endl;
    std::cout << '\n';

    // std::cout << "The total measured charge is: " << total_charge << " nC" << std::endl;
    // std::cout << '\n';

    std::cout << "The total neutron flux is: " << total_flux << " n cm^-2 s^-1" << std::endl;
    std::cout << "The uncertainty on the total flux is: " << total_flux_uncertainty << " n cm^-2 s^-1" << std::endl;
    std::cout << '\n';

    std::cout << "The average neutron energy is: " << avg_energy << " MeV" << std::endl;
    std::cout << "The uncertainty on the average energy is: " << avg_energy_uncertainty << " MeV" << std::endl;
    std::cout << '\n';

    // std::cout << "The neutron source strength is: " << source_strength << " n Gy^-1" << std::endl;
    // std::cout << '\n';


    //----------------------------------------------------------------------------------------------
    // Save spectrum to file
    //----------------------------------------------------------------------------------------------
    saveSpectrumAsRow(settings.path_output_spectra, num_bins, irradiation_conditions, spectrum, 
        spectrum_uncertainty, spectrum_uncertainty, energy_bins
    );
    std::cout << "Saved unfolded spectrum to " << settings.path_output_spectra << "\n";

    //----------------------------------------------------------------------------------------------
    // Generate report
    //----------------------------------------------------------------------------------------------
    if (settings.generate_report) {
        std::vector<double> measurements_report;
        if (settings.meas_units == "cps") {
            measurements_report = measurements;
            std::reverse(measurements_report.begin(),measurements_report.end());
        }
        else
            measurements_report = measurements_nc;

        UnfoldingReport myreport;

        if (settings.path_report.empty()) {
            settings.path_report = "output/report_" + irradiation_conditions + ".txt";
        }

        myreport.set_algorithm(settings.algorithm);
        myreport.set_path(settings.path_report);
        myreport.set_irradiation_conditions(irradiation_conditions);
        myreport.set_input_files(input_files);
        myreport.set_input_file_flags(input_file_flags);
        myreport.set_cutoff(settings.cutoff);
        myreport.set_error(settings.error);
        myreport.set_norm(settings.norm);
        myreport.set_f_factor(f_factor_report);
        myreport.set_num_measurements(num_measurements);
        myreport.set_num_bins(num_bins);
        myreport.set_num_poisson_samples(settings.num_poisson_samples);
        myreport.set_git_commit(GIT_COMMIT);
        myreport.set_measurements(measurements_report);
        myreport.set_dose_mu(dose_mu);
        myreport.set_doserate_mu(doserate_mu);
        myreport.set_duration(duration);
        myreport.set_meas_units(settings.meas_units);
        myreport.set_initial_spectrum(initial_spectrum);
        myreport.set_energy_bins(energy_bins);
        myreport.set_nns_response(nns_response);
        myreport.set_icrp_factors(icrp_factors);
        myreport.set_spectrum(spectrum);
        myreport.set_spectrum_uncertainty(spectrum_uncertainty);
        myreport.set_num_iterations(num_iterations);
        myreport.set_mlem_ratio(mlem_ratio);
        myreport.set_dose(ambient_dose_eq);
        myreport.set_s_dose(ambient_dose_eq_uncertainty);
        myreport.set_total_flux(total_flux);
        myreport.set_total_flux_uncertainty(total_flux_uncertainty);
        myreport.set_avg_energy(avg_energy);
        myreport.set_avg_energy_uncertainty(avg_energy_uncertainty);
        if (settings.algorithm == "mlemstop") {
            myreport.set_cps_crossover(settings.cps_crossover);
            myreport.set_j_threshold(j_threshold);
            myreport.set_j_final(j_factor);
        }
        myreport.prepare_report();

        std::cout << "Generated summary report: " << settings.path_report << "\n\n";
    }

    //----------------------------------------------------------------------------------------------
    // Plot the spectrum
    //----------------------------------------------------------------------------------------------
    if (settings.generate_figure) {
        std::cout << "Plotting spectrum: \n";
        plotSpectrum(figure_file_pre, figure_file_suf, irradiation_conditions, num_measurements, 
            num_bins, energy_bins, spectrum, spectrum_uncertainty
        );
        std::cout << "\n";
    }

    return 0;
}