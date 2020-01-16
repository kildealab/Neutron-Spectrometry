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
#include <chrono>
#include <random>
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
        "input/unfold_spectrum.cfg"
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

    double f_factor_report = settings.f_factor; // original value read in
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


    // if (settings.meas_units == "cps") {
    //     measurements = getMeasurements(settings);
    //     num_measurements = measurements.size();
    //     std::reverse(measurements.begin(),measurements.end());
    // }
    // // nC
    // else {
    //     measurements_nc = getMeasurements(settings);
    //     num_measurements = measurements_nc.size();
    //     // Convert to CPS
    //     for (int index=0; index < num_measurements; index++) {
    //         double measurement_cps = measurements_nc[num_measurements-index-1]*settings.norm/settings.f_factor/settings.duration;
    //         measurements.push_back(measurement_cps);
    //     }
    // }

    // Calculate average CPS value
    double avg_cps = 0;
    for (int i_meas = 0; i_meas < num_measurements; i_meas++) {
        avg_cps += measurements[i_meas];
    }
    avg_cps /= num_measurements;

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

    if (settings.num_meas_per_shell > 1) {
        std::cout << "The standard errors in CPS are:" << '\n'; // newline

        //Loop over the data matrix, display each value
        for (int i = 0; i < num_measurements; ++i)
        {
            std::cout << std_errors[i] << '\n';
        }
        std::cout << '\n';
    }

    // If doing MLEM-STOP: need ratio between "ideal" or "crossover" CPS value and the average
    // measured CPS. Call this scale_factor
    // double scale_factor = settings.cps_crossover/avg_cps;

    // Apply scaling factor to the measured values so that their average value equals the ideal CPS
    // Then MLEM-STOP can be done such that terminate when J=1 (rather than J=scale_factor)
    // for (int i_meas = 0; i_meas < num_measurements; i_meas++) {
    //     measurements[i_meas] *= scale_factor;
    // }

    //----------------------------------------------------------------------------------------------
    // Generate the energy bins matrix:
    //  - size = # of energy bins
    // Input the energies from energy bins file
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

    // MLEM-STOP specific parameters, initialized here for use later
    double j_factor = 0;
    double j_threshold = 0;

    // Unfold spectrum according to user-specified algorithm
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
        throw std::logic_error("Unrecognized unfolding algorithm: " + settings.algorithm);
        // std::cout << "No unfolding algorithm found for: " + settings.algorithm + '\n';
    }

    // Need to "unscale" spectrum back to true values in order to calculate quantities of interest
    // for (int i_bin = 0; i_bin < num_bins; i_bin++) {
    //     spectrum[i_bin] /= scale_factor;
    // }

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
        // std::cout << mlem_estimate[i_meas]/scale_factor << '\n';
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
    if (settings.algorithm == "mlemstop") {
        std::cout << "J factor: " << j_factor << "\n";
        std::cout << "J threshold: " << j_threshold << "\n";
    }

    //----------------------------------------------------------------------------------------------
    // Calculate quantities of interest (e.g. dose & its uncertainty)
    //----------------------------------------------------------------------------------------------
    double ambient_dose_eq = calculateDose(num_bins, spectrum, icrp_factors);
    // double total_charge = calculateTotalCharge(num_measurements,measurements_nc);
    double total_flux = calculateTotalFlux(num_bins,spectrum);
    double avg_energy = calculateAverageEnergy(num_bins,spectrum,energy_bins);
    // double source_strength = calculateSourceStrength(num_bins,spectrum,duration,dose_mu);

    // Need to scale spectrum again in order to do uncertainty calculations
    // for (int i_bin = 0; i_bin < num_bins; i_bin++) {
    //     spectrum[i_bin] *= scale_factor;
    // }

    //----------------------------------------------------------------------------------------------
    // Determine the uncertainty in the unfolded spectrum using one of the available methods
    //----------------------------------------------------------------------------------------------
    std::vector<double> spectrum_uncertainty;
    std::vector<double> spectrum_uncertainty_lower;
    std::vector<double> spectrum_uncertainty_upper;
    double ambient_dose_eq_uncertainty_upper = 0;
    double ambient_dose_eq_uncertainty_lower = 0;

    // MLEM-STOP specific parameters, initialized here for use later
    UncertaintyManagerJ j_manager_low(j_threshold,1+settings.sigma_j);
    UncertaintyManagerJ j_manager_high(j_threshold,1-settings.sigma_j);
    // The # of sampled measurement sets that are discarded b/c don't converge with MLEM-STOP
    int num_toss = 0;

    // This approach generates a series of sampled measurements (using original measurements as the
    // means). Unfolding is performed for each of these spectra. The uncertainty in the unfolded
    // spectrum is taken to be the Root-Mean-Square-Deviation between the unfolded spectrum and each
    // sampled spectrum. The number of samples is set by the user via num_uncertainty_samples
    if (settings.uncertainty_type == "poisson" || settings.uncertainty_type == "gaussian") {
        std::vector<std::vector<double>> sampled_spectra; // dimensions: num_uncertainty_samples x num_bins
        std::vector<double> sampled_dose; // dimension: num_uncertainty_samples

        for (int i_samp = 0; i_samp < settings.num_uncertainty_samples; i_samp++) {
            std::vector<double> sampled_measurements; // dimension: num_measurements
            std::vector<double> sampled_mlem_ratio; // dimension: num_measurements
            std::vector<double> sampled_mlem_correction; // dimension: num_measurements
            std::vector<double> sampled_mlem_estimate; // dimension: num_measurements
            std::vector<double> sampled_spectrum = initial_spectrum; // dimension: num_bins

            bool toss_sample = false; // Track if sample was tossed

            // If doing Poisson-sampling to generate pseudo-measurement set:
            if (settings.uncertainty_type == "poisson") {
                for (int i_meas = 0; i_meas < num_measurements; i_meas++) {
                    double sampled_value = 0;
                    for (int i_samp =0; i_samp < settings.num_meas_per_shell; i_samp++) {
                        sampled_value += poisson(measurements[i_meas]);
                    }
                    sampled_value /= settings.num_meas_per_shell;
                    sampled_measurements.push_back(sampled_value);
                }
            }
            // If doing Gaussian-sampling to generate pseudo-measurement set
            else if (settings.uncertainty_type == "gaussian") {
                for (int i_meas = 0; i_meas < num_measurements; i_meas++) {
                    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                    std::default_random_engine generator (seed);
                    std::normal_distribution<double> distribution(measurements[i_meas],std_errors[i_meas]);

                    double new_sample = distribution(generator);
                    sampled_measurements.push_back(new_sample);
                }
            }

            // Do unfolding on the initial spectrum & sampled measurement values
            if (settings.algorithm == "mlem") {
                runMLEM(settings.cutoff, settings.error, num_measurements, num_bins, sampled_measurements, 
                    sampled_spectrum, nns_response, normalized_response, sampled_mlem_ratio, 
                    sampled_mlem_correction, sampled_mlem_estimate
                );
            }
            // MLEM-STOP requires special handling of unfolding the sampled measurement sets.
            // Despite best efforts, sometimes MLEM-STOP will never converge for some samples. 
            // Current best approach is to discard those samples. A record is kept of the number
            // of samples discarded and reported to the user so they may interpret the final
            // uncertainty accordingly.
            else if (settings.algorithm == "mlemstop") {
                // Calculate unique J threshold for the current sample
                double sampled_j_threshold = determineJThreshold(num_measurements,sampled_measurements,settings.cps_crossover);
                double sampled_j_factor = 0;

                // Try running MLEM-STOP for the current sample but catch an error if never
                // converges. Increment num_toss if so, and decrement the sample index (i.e scrap
                // the current sample and retry)
                try {
                    runMLEMSTOP(settings.cutoff, num_measurements, num_bins, sampled_measurements,
                        sampled_spectrum, nns_response, normalized_response, sampled_mlem_ratio, sampled_mlem_correction, 
                        sampled_mlem_estimate, sampled_j_threshold, sampled_j_factor
                    );
                }
                catch (std::logic_error e) {
                    i_samp = i_samp - 1;
                    num_toss = num_toss + 1;
                    toss_sample = true;
                }
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
                throw std::logic_error("Unrecognized unfolding algorithm: " + settings.algorithm);
            }

            // Note the sampled CPS values were based on the scaled CPS, so need to scale back the
            // sampled spectra to correct magnitude 
            // for (int i_bin = 0; i_bin < num_bins; i_bin++) {
            //     sampled_spectrum[i_bin] /= scale_factor;
            // }

            if (!toss_sample) {
                sampled_spectra.push_back(sampled_spectrum); // add to growing array of sampled spectra
            }

            // Calculate the ambient dose equivalent associated with the sampled spectrum
            double sdose = calculateDose(num_bins, sampled_spectrum, icrp_factors);

            sampled_dose.push_back(sdose);
        }

        // Finally, "unscale" spectrum back to true values for remaining calculations & logging
        // for (int i_bin = 0; i_bin < num_bins; i_bin++) {
        //     spectrum[i_bin] /= scale_factor;
        // }

        // Calculate the spectrum uncertainty (same upper & lower)
        calculateRMSD_vector(settings.num_uncertainty_samples, spectrum, sampled_spectra, spectrum_uncertainty_lower);
        spectrum_uncertainty_upper = spectrum_uncertainty_lower;

        // Calculate the dose uncertainty (same upper & lower)
        ambient_dose_eq_uncertainty_upper = calculateRMSD(settings.num_uncertainty_samples, ambient_dose_eq, sampled_dose);
        ambient_dose_eq_uncertainty_lower = ambient_dose_eq_uncertainty_upper;

        // If want to print the number of sample sets kept vs tossed:
        // std::cout << "Number of sampled measurement sets kept: " << settings.num_uncertainty_samples << "\n";
        // std::cout << "Number of sampled measurement sets tossed: " << num_toss << "\n";
    }

    // if (settings.uncertainty_type == "gaussian") {
    //     std::vector<std::vector<double>> sampled_spectra; // dimensions: num_uncertainty_samples x num_bins
    //     std::vector<double> sampled_dose; // dimension: num_uncertainty_samples


    // }

    // This approach uses a known uncertainty around J=1, specified by sigma_j. An "upper" and "lower"
    // spectrum estimate are generated using the user-defined sigma_j parameter. A unique j_threshold
    // is determined for both upper and lower, and MLEM-STOP is performed using both thresholds. The
    // upper uncertainty is then the difference betwen the "upper" spectrum and the MLEM-STOP estimated 
    // spectrum. Similarly for the lower uncertainty. This is all handled in the UncertaintyManagerJ
    // class.
    else if (settings.uncertainty_type == "j_bounds") {
        j_manager_low.determineSpectrumUncertainty(spectrum,settings.cutoff,num_measurements,
            num_bins,measurements,nns_response,normalized_response,initial_spectrum
        );
        spectrum_uncertainty_lower = j_manager_low.spectrum_uncertainty;

        j_manager_high.determineSpectrumUncertainty(spectrum,settings.cutoff,num_measurements,
            num_bins,measurements,nns_response,normalized_response,initial_spectrum
        );
        spectrum_uncertainty_upper = j_manager_high.spectrum_uncertainty;

        j_manager_low.determineDoseUncertainty(ambient_dose_eq,spectrum,num_bins,icrp_factors);
        ambient_dose_eq_uncertainty_lower = j_manager_low.dose_uncertainty;

        j_manager_high.determineDoseUncertainty(ambient_dose_eq,spectrum,num_bins,icrp_factors);
        ambient_dose_eq_uncertainty_upper = j_manager_high.dose_uncertainty;
    }
    else {
        throw std::logic_error("Unrecognized uncertainty type: " + settings.uncertainty_type);
    }

    //----------------------------------------------------------------------------------------------
    // Calculate uncertainties in quantities of interest
    //----------------------------------------------------------------------------------------------
    double total_flux_uncertainty_upper = calculateSumUncertainty(num_bins,spectrum_uncertainty_upper);
    double total_flux_uncertainty_lower = calculateSumUncertainty(num_bins,spectrum_uncertainty_lower);

    double avg_energy_uncertainty_upper = calculateEnergyUncertainty(num_bins,energy_bins,spectrum,
        spectrum_uncertainty_upper,total_flux,total_flux_uncertainty_upper
    );
    double avg_energy_uncertainty_lower = calculateEnergyUncertainty(num_bins,energy_bins,spectrum,
        spectrum_uncertainty_lower,total_flux,total_flux_uncertainty_lower
    );

    //----------------------------------------------------------------------------------------------
    // Display calculated quantities
    //----------------------------------------------------------------------------------------------
    std::cout << '\n';
    std::cout << "The equivalent dose is: " << ambient_dose_eq << " mSv/h" << std::endl;
    std::cout << "Upper uncertainty: " << ambient_dose_eq_uncertainty_upper << " mSv/h" << std::endl;
    std::cout << "Lower uncertainty: " << ambient_dose_eq_uncertainty_lower << " mSv/h" << std::endl;
    std::cout << '\n';

    // std::cout << "The total measured charge is: " << total_charge << " nC" << std::endl;
    // std::cout << '\n';

    std::cout << "The total neutron flux is: " << total_flux << " n cm^-2 s^-1" << std::endl;
    std::cout << "Upper uncertainty: " << total_flux_uncertainty_upper << " n cm^-2 s^-1" << std::endl;
    std::cout << "Lower uncertainty: " << total_flux_uncertainty_lower << " n cm^-2 s^-1" << std::endl;
    std::cout << '\n';

    std::cout << "The average neutron energy is: " << avg_energy << " MeV" << std::endl;
    std::cout << "Upper uncertainty: " << avg_energy_uncertainty_upper << " MeV" << std::endl;
    std::cout << "Lower uncertainty: " << avg_energy_uncertainty_lower << " MeV" << std::endl;

    std::cout << '\n';

    // std::cout << "The neutron source strength is: " << source_strength << " n Gy^-1" << std::endl;
    // std::cout << '\n';


    //----------------------------------------------------------------------------------------------
    // Save spectrum to file
    //----------------------------------------------------------------------------------------------
    saveSpectrumAsRow(settings.path_output_spectra, num_bins, settings.irradiation_conditions, spectrum, 
        spectrum_uncertainty_upper, spectrum_uncertainty_lower, energy_bins
    );
    std::cout << "Saved unfolded spectrum to " << settings.path_output_spectra << "\n";

    //----------------------------------------------------------------------------------------------
    // Generate report
    //----------------------------------------------------------------------------------------------
    if (settings.generate_report) {
        // std::vector<double> measurements_report;
        // if (settings.meas_units == "cps") {
        //     measurements_report = measurements;
        //     std::reverse(measurements_report.begin(),measurements_report.end());
        // }
        // else
        //     measurements_report = measurements_nc;

        UnfoldingReport myreport;

        if (settings.path_report.empty()) {
            settings.path_report = "output/report_" + settings.irradiation_conditions + ".txt";
        }

        myreport.set_algorithm(settings.algorithm);
        myreport.set_path(settings.path_report);
        myreport.set_irradiation_conditions(settings.irradiation_conditions);
        myreport.set_input_files(input_files);
        myreport.set_input_file_flags(input_file_flags);
        myreport.set_cutoff(settings.cutoff);
        myreport.set_error(settings.error);
        myreport.set_norm(settings.norm);
        myreport.set_f_factor(f_factor_report);
        myreport.set_num_measurements(num_measurements);
        myreport.set_uncertainty_type(settings.uncertainty_type);
        myreport.set_num_bins(num_bins);
        myreport.set_num_uncertainty_samples(settings.num_uncertainty_samples);
        myreport.set_git_commit(GIT_COMMIT);
        myreport.set_measurements(measurements);
        myreport.set_measurements_nc(measurements_nc);
        myreport.set_dose_mu(settings.dose_mu);
        myreport.set_doserate_mu(settings.doserate_mu);
        myreport.set_duration(settings.duration);
        myreport.set_meas_units(settings.meas_units);
        myreport.set_initial_spectrum(initial_spectrum);
        myreport.set_energy_bins(energy_bins);
        myreport.set_nns_response(nns_response);
        myreport.set_icrp_factors(icrp_factors);
        myreport.set_spectrum(spectrum);
        myreport.set_spectrum_uncertainty_upper(spectrum_uncertainty_upper);
        myreport.set_spectrum_uncertainty_lower(spectrum_uncertainty_lower);
        myreport.set_num_iterations(num_iterations);
        myreport.set_mlem_ratio(mlem_ratio);
        myreport.set_dose(ambient_dose_eq);
        myreport.set_dose_uncertainty_upper(ambient_dose_eq_uncertainty_upper);
        myreport.set_dose_uncertainty_lower(ambient_dose_eq_uncertainty_lower);
        myreport.set_total_flux(total_flux);
        myreport.set_total_flux_uncertainty_upper(total_flux_uncertainty_upper);
        myreport.set_total_flux_uncertainty_lower(total_flux_uncertainty_lower);
        myreport.set_avg_energy(avg_energy);
        myreport.set_avg_energy_uncertainty_upper(avg_energy_uncertainty_upper);
        myreport.set_avg_energy_uncertainty_lower(avg_energy_uncertainty_lower);
        if (settings.algorithm == "mlemstop") {
            myreport.set_cps_crossover(settings.cps_crossover);
            myreport.set_j_threshold(j_threshold);
            myreport.set_j_final(j_factor);
            myreport.set_j_manager_low(j_manager_low);
            myreport.set_j_manager_high(j_manager_high);
            myreport.set_num_toss(num_toss);
        }
        myreport.prepare_report();

        std::cout << "Generated summary report: " << settings.path_report << "\n\n";
    }

    //----------------------------------------------------------------------------------------------
    // Plot the spectrum
    //----------------------------------------------------------------------------------------------
    if (settings.generate_figure) {
        std::cout << "Plotting spectrum: \n";
        if (settings.path_figure.empty()) {
            settings.path_figure = "output/figure_" + settings.irradiation_conditions + ".png";
        }
        plotSpectrum(settings.path_figure, settings.irradiation_conditions, num_measurements, 
            num_bins, energy_bins, spectrum, spectrum_uncertainty_upper, spectrum_uncertainty_lower
        );
        std::cout << "\n";
    }

    return 0;
}