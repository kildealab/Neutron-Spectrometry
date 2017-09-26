//**************************************************************************************************
// This code calls the files bins.csv, respmat.csv and ini.csv and uses the MLEM algorithm to create an output matrix
// of dimensions 52*1 from a data matrix of dimensions 8*1. In addition, Poisson distribution is used to calculate a
// pseudo measured matrix of dimensions 8*1, then MLEM algorithm is used to calculate a pseudo output matrix of dimensions
// 52*1; this procedure is done 1000 times. After this, an ini_poisson matrix is created using the 1000 pseudo output matrices
// resulting in a matrix of dimensions 52*1000.
// Also, this code calculates the standard deviation and the variance by using the elements from the two matrices
// (52*1) and (52*1000).
// In addition, this code uses the ROOT libraries to plot an histogram of the neutron fluence rate (the output matrix of
// dimensions 52*1) in function of the energy bins (The result in an image.png showing the graph).
//**************************************************************************************************
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <random>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <algorithm>

// Root
#include "TAxis.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TPad.h"
#include "TROOT.h"
#include "TColor.h"
#include "TFrame.h"
#include "TVirtualPad.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TText.h"

#include <vector>

// Prototypes
double poisson(double lambda);
bool is_empty(std::ifstream& pFile);

int setfile(std::vector<std::string> &arg_vector, std::string directory, std::string arg_string, std::string default_filename, std::string &filename);
void checkUnknownParameters(std::vector<std::string> &arg_vector, std::vector<std::string> input_file_flags);
int setSettings(std::string config_file, int &cutoff, double &norm, double &error, double &f_factor, int &num_poisson_samples);
std::vector<double> getMeasurements(std::string input_file, std::string &irradiation_conditions, double &dose_mu, double &doserate_mu, int &t);
int saveDose(std::string dose_file, std::string irradiation_conditions, double dose, double s_dose);
int saveSpectrum(std::string spectrum_file, std::string irradiation_conditions, std::vector<double>& spectrum, std::vector<double>& spectrum_uncertainty, std::vector<double>& energy_bins);
int prepareReport(std::string report_file, std::string irradiation_conditions, std::vector<std::string> &input_files, std::vector<std::string> &input_file_flags, int cutoff, double error, double norm, double f_factor, int num_measurements, int num_bins, int num_poisson_samples, std::vector<double>& measurements_nc, double dose_mu, double doserate_mu, int duration, std::vector<double>& energy_bins, std::vector<double>& initial_spectrum, std::vector<std::vector<double>>& nns_response, int num_iterations, std::vector<double>& mlem_ratio, double dose, double s_dose, std::vector<double>& spectrum, std::vector<double>& spectrum_uncertainty, std::vector<double>& icrp_factors);
int readInputFile1D(std::string file_name, std::vector<double>& input_vector);
int readInputFile2D(std::string file_name, std::vector<std::vector<double>>& input_vector);
int checkDimensions(int reference_size, std::string reference_string, int test_size, std::string test_string);
int runMLEM(int cutoff, double error, int num_measurements, int num_bins, std::vector<double> &measurements, std::vector<double> &spectrum, std::vector<std::vector<double>>& nns_response, std::vector<double> &mlem_ratio);
double calculateDose(int num_bins, std::vector<double> &spectrum, std::vector<double> &icrp_factors);
int calculateRMSD_vector(int num_samples, std::vector<double> &true_vector, std::vector<std::vector<double>> &sampled_vectors, std::vector<double> &rms_differences);
double calculateRMSD(int num_samples, double true_value, std::vector<double> &sample_vector);
int plotSpectrum(std::string figure_file_pre, std::string figure_file_suf, std::string irradiation_conditions, int num_measurements, int num_bins, std::vector<double> &energy_bins, std::vector<double> &spectrum, std::vector<double> &spectrum_uncertainty);

// Constants
std::string DOSE_HEADERS[] = {
    "Irradiation Conditions",
    "Dose Rate (mSv/hr)",
    "RMS Error (mSv/hr)"
};

std::string UNCERTAINTY_SUFFIX = "_ERROR";

// mt19937 is a random number generator class based on the Mersenne Twister algorithm
// Create an mt19937 object, called mrand, that is seeded with the current time in seconds
std::mt19937 mrand(std::time(0));

//==================================================================================================
//==================================================================================================
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
        double measurement_cps = measurements_nc[num_measurements-index-1]*norm/f_factor/duration*(dose_mu/doserate_mu);
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
    prepareReport(report_file, irradiation_conditions, input_files, input_file_flags, cutoff, error, norm, f_factor_report, num_measurements, num_bins, num_poisson_samples, measurements_nc, dose_mu, doserate_mu, duration, energy_bins, initial_spectrum, nns_response, num_iterations, mlem_ratio, ambient_dose_eq, ambient_dose_eq_uncertainty, spectrum, spectrum_uncertainty, icrp_factors);
    std::cout << "Generated summary report: " << report_file << "\n\n";

    //----------------------------------------------------------------------------------------------
    // Plot the spectrum
    //----------------------------------------------------------------------------------------------
    plotSpectrum(figure_file_pre, figure_file_suf, irradiation_conditions, num_measurements, num_bins, energy_bins, spectrum, spectrum_uncertainty);

    return 0;
}

//**************************************************************************************************
// Helper functions
//**************************************************************************************************

//==================================================================================================
// Return a poisson_distribution object, d that can be sampled from
//  - accept a parameter 'lamda', that is used as the mean & std. dev of the poisson distribution
//  - use the 'mrand' generator as the random number generator for the distribution
//==================================================================================================
double poisson(double lambda)
{
    std::poisson_distribution<int> d(lambda); // initialization

    return d(mrand); // sample
}

//==================================================================================================
// Determine if a file is empty
//==================================================================================================
bool is_empty(std::ifstream& pFile)
{
    return pFile.peek() == std::ifstream::traits_type::eof();
}

//==================================================================================================
// Check arguments passed to the main function for a matching flag, and assign value for the input
// filename
// Args:
//  - arg_vector: Vector containing all arguments passed to the main function (at command line)
//  - directory: Directory (path) at which the input file is located (e.g. input/)
//  - arg_string: The argument string that indicates what the file is (e.g. --measurements)
//  - default_file: The default filename to be used for the arg_string provided
//  - filename: The variable (passed by reference) that will be assigned a value, representing the
//      the filename to be used
//==================================================================================================
int setfile(std::vector<std::string> &arg_vector, std::string directory, std::string arg_string, std::string default_filename, std::string &filename) {
    // Place an iterator at an index of the vector, where the matching arg_string was found
    std::vector<std::string>::iterator iter_args = std::find(arg_vector.begin(), arg_vector.end(), arg_string);

    // If a matching argument was found, ensure that a subsequent argument (i.e. filename to use) 
    // was provided. Apply it if so. If not, throw an error
    if( iter_args != arg_vector.end()) {
        iter_args += 1;
        if (iter_args != arg_vector.end()) {
            filename = directory + *iter_args;
        }
        else {
            // throw error saying no file provided for *iter_measurements -= 1
            iter_args -= 1;
            std::cout << "Error: no file provided for argument: " << *iter_args << "\n";
        }
    }
    // If no match was found for the target arg_string within arg_vector, use the default filename
    // i.e. The user did not specify which file to use
    else {
        filename = directory + default_filename;
    }

    return 1;
}

//==================================================================================================
// Identify any options that were provided to the main function that are not supported and notify
// the user.
//==================================================================================================
void checkUnknownParameters(std::vector<std::string> &arg_vector, std::vector<std::string> input_file_flags) {
    for (int i=0; i<arg_vector.size(); i++) {
        std::vector<std::string>::iterator iter_args = std::find(input_file_flags.begin(), input_file_flags.end(), arg_vector[i]);
        if(iter_args == input_file_flags.end()) {
            std::cout << "Warning: Ignored unknown argument " << arg_vector[i] << "\n";
        }
        else {
            i += 1;
        }
    }
   
}

//==================================================================================================
// Retrieve settings from a configuration file 'config_file' and save values in relevant variables
//==================================================================================================
int setSettings(std::string config_file, int &cutoff, double &norm, double &error, double &f_factor, int &num_poisson_samples) {
    std::ifstream cfile(config_file);
    std::string line;
    std::vector<double> settings;

    // If file is able to be read
    if (cfile.is_open())
    {
        // loop through each line in the file, extract the value for each setting into 'token'
        while ( getline (cfile,line) )
        {
            // settings format: setting_name=value
            std::string delimiter = "=";
            std::string token = line.substr(line.find(delimiter)+1); // substring from '=' to end of string
            // std::cout << token << '\n';
            settings.push_back(atof(token.c_str())); // convert str to double, insert into vector
        }
        cfile.close();
    }
    // Problem opening the file
    else {
        return false;
    }

    // Assign the settings values to function parameters passed by reference
    cutoff = (int)settings[0];
    norm = settings[1];
    error = settings[2];
    f_factor = settings[3];
    num_poisson_samples = settings[4];
    return true;
}

//==================================================================================================
// Read measurement data from file
// data_vector is ordered from index:0 storing the value for 7 moderators to index:7 storing the value
// for 0 moderators
//==================================================================================================
std::vector<double> getMeasurements(std::string input_file, std::string &irradiation_conditions, double &dose_mu, double &doserate_mu, int &t) {
    std::ifstream ifile(input_file);
    if (!ifile.is_open()) {
        //throw error
        std::cout << "Unable to open input file:" + input_file + '\n';
    }

    // Load header information from 'ifile'
    getline(ifile,irradiation_conditions);
    // removes: carriage return '\r' from the string (which causes weird string overwriting)
    irradiation_conditions.erase( std::remove(irradiation_conditions.begin(), irradiation_conditions.end(), '\r'), irradiation_conditions.end() );

    // Extract dose & measurement duration
    std::string dose_string;
    getline(ifile,dose_string);
    dose_mu = atoi(dose_string.c_str());

    std::string doserate_string;
    getline(ifile,doserate_string);
    doserate_mu = atoi(doserate_string.c_str());

    std::string t_string;
    getline(ifile,t_string);
    t = atoi(t_string.c_str());

    // Loop through file, get measurement data
    std::string line;
    std::vector<double> data_vector;
    while (getline(ifile,line)) {
        std::istringstream line_stream(line);
        std::string stoken; // store individual values between delimiters on a line

        // Loop through each line, delimiting at commas
        while (getline(line_stream, stoken, ',')) {
            data_vector.push_back(atof(stoken.c_str())); // add data to the vector
        }
    }

    // print elements of data vector
    // int vector_size = data_vector.size();
    // for (int i=0; i<vector_size; i++) {
    //     std::cout << data_vector[i] << '\n';
    // }
    ifile.close();
    std::cout << "Data successfully retrieved from " + input_file + '\n';
    return data_vector;
}

//==================================================================================================
// Save calculated dose (and its error) to file
//==================================================================================================
int saveDose(std::string dose_file, std::string irradiation_conditions, double dose, double s_dose) {
    // determine if file exists
    std::ifstream checkfile(dose_file);
    bool file_empty = is_empty(checkfile);
    checkfile.close();

    // If file does not exist, create it. If file exists, open it
    std::ofstream dfile;
    dfile.open(dose_file, std::ios_base::app);

    // Add header line to start of file, if file was empty
    if (file_empty) {
        std::ostringstream header_stream;
        header_stream << DOSE_HEADERS[0] << "," << DOSE_HEADERS[1] << "," << DOSE_HEADERS[2] << "\n";
        std::string header = header_stream.str();
        dfile << header; 
    }

    // Append new line of data
    std::ostringstream new_data_stream;
    new_data_stream << irradiation_conditions << "," << dose << "," << s_dose << "\n";
    std::string new_data = new_data_stream.str();
    dfile << new_data; 

    dfile.close();
    return 1;
}

//==================================================================================================
// Save calculated spectrum (and error spectrum) to file
//==================================================================================================
int saveSpectrum(std::string spectrum_file, std::string irradiation_conditions, std::vector<double>& spectrum, std::vector<double> &spectrum_uncertainty, std::vector<double>& energy_bins) {
    // determine if file exists
    std::ifstream sfile(spectrum_file);
    bool file_empty = is_empty(sfile);
    bool file_exists = sfile.good();

    std::vector<std::string> sfile_lines;

    // If the file exists and is not empty: append 2 new columns to existing rows (strings)
    // First new column: the spectrum, Second new column: the uncertainty on the spectrum
    if (file_exists && !file_empty) {
        // Retrieve and update the header
        int index = 0;
        std::string header;
        getline(sfile, header);
        header = header + "," + irradiation_conditions + "," + irradiation_conditions + UNCERTAINTY_SUFFIX + "\n";
        sfile_lines.push_back(header);

        // Loop through existing lines in file (including header)
        // For each, string-concatenate the new data
        // Store new lines in vector 'sfile_lines'
        std::string line;
        while (getline(sfile,line)) {
            //Note: using stream b/c have to convert doubles to strings
            std::ostringstream line_stream;
            line_stream << line << "," << spectrum[index] << "," << spectrum_uncertainty[index] << "\n";
            line = line_stream.str();
            sfile_lines.push_back(line);
            index++;
        }
    }
    // If the file does not exist or was empty: create each line for the file
    // Column 1: energy bins, Column 2: the spectrum, Column 3: the uncertainty on spectrum
    else {
        // Prepare header
        std::string header;
        header = "Energy (MeV)," + irradiation_conditions + "," + irradiation_conditions + UNCERTAINTY_SUFFIX + "\n";
        sfile_lines.push_back(header);

        // Prepare contents
        std::string line;
        for (int index=0; index<spectrum.size(); index++) {
            //Note: using stream b/c have to convert doubles to strings
            std::ostringstream line_stream;
            line_stream << energy_bins[index] << "," << spectrum[index] << "," << spectrum_uncertainty[index] << "\n";
            line = line_stream.str();
            sfile_lines.push_back(line);
        }
    }
    sfile.close();

    // Rewrite file using lines stored in vector 'sfile_lines'
    std::ofstream nfile(spectrum_file);
    int vector_size = sfile_lines.size();
    for (int i=0; i<vector_size; i++) {
        nfile << sfile_lines[i];
    }

    return 1;
}

//==================================================================================================
// Generate a textfile report of pertinent information from execution of this program. The contents
// of this function are separated by headers indicating the type of information printed to the
// report in the the corresponding section.
//==================================================================================================
int prepareReport(std::string report_file, std::string irradiation_conditions, std::vector<std::string> &input_files, std::vector<std::string> &input_file_flags, int cutoff, double error, double norm, double f_factor, int num_measurements, int num_bins, int num_poisson_samples, std::vector<double>& measurements_nc, double dose_mu, double doserate_mu, int duration, std::vector<double>& energy_bins, std::vector<double>& initial_spectrum, std::vector<std::vector<double>>& nns_response, int num_iterations, std::vector<double>& mlem_ratio, double dose, double s_dose, std::vector<double>& spectrum, std::vector<double>& spectrum_uncertainty, std::vector<double>& icrp_factors) {
    std::string HEADER_DIVIDE = "************************************************************************************************************************\n";
    std::string SECTION_DIVIDE = "\n========================================================================================================================\n\n";
    std::string COLSTRING = "--------------------";
    int sw = 30; // settings column width
    int cw = 20; // data column width
    int rw = 9; // NNS response column width

    std::ofstream rfile(report_file);

    //----------------------------------------------------------------------------------------------
    // Header
    //----------------------------------------------------------------------------------------------
    rfile << HEADER_DIVIDE;
    rfile << "Neutron Spectrometry Report\n\n";
    rfile << std::left << std::setw(sw) << "Irradiation Specs: " << irradiation_conditions << "\n";
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    rfile << std::left << std::setw(sw) << "Date report was generated: " << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << "\n";
    rfile << "Input arguments (files) used:\n";
    for (int i=0; i<input_files.size(); i++) {
        std::string tempstring = "    " + input_file_flags[i];
        rfile << std::left <<std::setw(sw) << tempstring << input_files[i] << "\n";
    }
    rfile << HEADER_DIVIDE << "\n";

    //----------------------------------------------------------------------------------------------
    // Settings
    //----------------------------------------------------------------------------------------------
    rfile << "Settings\n\n";
    rfile << std::left << std::setw(sw) << "MLEM max # of iterations:" << cutoff << "\n";
    rfile << std::left << std::setw(sw) << "MLEM target ratio:" << error << "\n";
    rfile << std::left << std::setw(sw) << "NNS normalization factor:" << norm << "\n";
    rfile << std::left << std::setw(sw) << "NNS calibration factor:" << f_factor << " fA/cps\n";
    rfile << std::left << std::setw(sw) << "Number of poisson samples:" << num_poisson_samples << "\n";
    rfile << SECTION_DIVIDE;

    //----------------------------------------------------------------------------------------------
    // Measurement
    //----------------------------------------------------------------------------------------------
    rfile << "Measurement\n\n";
    rfile << std::left << std::setw(sw) << "Delivered dose:" << dose_mu << " MU\n";
    rfile << std::left << std::setw(sw) << "Delivered doserate:" << doserate_mu << " MU/min\n";
    rfile << std::left << std::setw(sw) << "Measurement duration:" << duration << " s\n\n";
    // rfile << "Measured Data (measurement duration: " << duration << "s)\n\n";
    rfile << std::left << std::setw(cw) << "# of moderators" << "Charge (nC)\n";
    rfile << std::left << std::setw(cw) << COLSTRING << COLSTRING << "\n";
    for (int i=0; i<num_measurements; i++) {
        rfile << std::left << std::setw(cw) << i << measurements_nc[i] << "\n";
    }
    rfile << SECTION_DIVIDE;

    //----------------------------------------------------------------------------------------------
    // Inputs
    //----------------------------------------------------------------------------------------------
    rfile << "Inputs (Number of energy bins: " << num_bins << ")\n\n";
    rfile << std::left << std::setw(cw) << "Energy bins" << std::setw(cw) << "Input spectrum" << "| NNS Response by # of moderators (cm^2)\n";
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

    //----------------------------------------------------------------------------------------------
    // MLEM Processing
    //----------------------------------------------------------------------------------------------
    rfile << "MLEM information\n\n";
    rfile << std::left << std::setw(sw) << "# of iterations: " << num_iterations << "/" << cutoff << "\n\n";
    rfile << "Final MLEM ratio = measured charge / estimated charge:\n";
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

    //----------------------------------------------------------------------------------------------
    // Results
    //----------------------------------------------------------------------------------------------
    rfile << "Results\n\n";
    rfile << std::left << std::setw(sw) << "Ambient dose equivalent:" << dose << " mSv/hr\n";
    rfile << std::left << std::setw(sw) << "Uncertainty:" << s_dose << " mSv\n\n";
    rfile << std::left << std::setw(cw) << "Energy bins" << std::setw(cw) << "Unfolded spectrum" << std::setw(cw) << "Uncertainty" << std::setw(cw) << "| ICRP H factor" << "Ambient Dose Equiv.\n";
    rfile << std::left << std::setw(cw) << "(MeV)" << std::setw(cw) << "(n cm^-2 s^-1)" << std::setw(cw) << "(n cm^-2 s^-1)" << std::setw(cw) << "| (pSv/cm^2)" << "(mSv/hr)\n";;
    rfile << std::left << std::setw(cw) << COLSTRING << std::setw(cw) << COLSTRING << std::setw(cw) << COLSTRING << std::setw(cw) << COLSTRING << COLSTRING << "\n";
    for (int i=0; i<num_bins; i++) {
        std::ostringstream icrp_string;
        icrp_string << "| " <<icrp_factors[i];
        double subdose = spectrum[i]*icrp_factors[i]*3600*(1e-9);
        rfile << std::left << std::setw(cw) << energy_bins[i] << std::setw(cw) << spectrum[i] << std::setw(cw) << spectrum_uncertainty[i] << std::setw(26) << icrp_string.str() << subdose << "\n";
    }

    rfile.close();
    return 1;
}

//==================================================================================================
// Read a 1D CSV file (i.e. one value per row) and store the data in a vector
//==================================================================================================
int readInputFile1D(std::string file_name, std::vector<double>& input_vector) {
    std::ifstream ifile(file_name);
    std::string iline;

    if (!ifile.is_open()) {
        //throw error
        throw std::logic_error("Unable to open input file: " + file_name);
    }
    while (getline(ifile,iline)) {
        std::istringstream line_stream(iline);
        std::string stoken; // store individual values between delimiters on a line

        // Delimit line at trailing comma
        getline(line_stream, stoken, ',');
        input_vector.push_back(atof(stoken.c_str()));
    }
    return 1;
}

//==================================================================================================
// Read a CSV file and store the data in a 2D vector
//==================================================================================================
int readInputFile2D(std::string file_name, std::vector<std::vector<double>>& input_vector) {
    std::ifstream ifile(file_name);
    std::string iline;

    if (!ifile.is_open()) {
        //throw error
        throw std::logic_error("Unable to open input file: " + file_name);
    }

    // Loop through each line in the file
    while (getline(ifile,iline)) {
        std::vector<double> new_column;

        std::istringstream line_stream(iline);
        std::string stoken; // store individual values between delimiters on a line

        // Loop through the line, delimiting at commas
        while (getline(line_stream, stoken, ',')) {
            new_column.push_back(atof(stoken.c_str())); // add data to the vector
        }
        input_vector.push_back(new_column);
    }
    return 1;
}

//==================================================================================================
// Compare the size of two vectors (indicated by reference_size and test_size). If the size of test
// does not match the size of reference, throw an error with a message constructed using text passed
// as function parameters (reference_string and test_string)
//==================================================================================================
int checkDimensions(int reference_size, std::string reference_string, int test_size, std::string test_string) {
    std::ostringstream error_message;
    if (reference_size != test_size) {
        error_message << "File dimensions mismatch: " << test_string << " (" << test_size << ") does not match " << reference_string << " (" << reference_size << ")";
        throw std::logic_error(error_message.str());   
    }
    return 1;
}

//==================================================================================================
// Accept a series of measurements and an estimated input spectrum and perform the MLEM algorithm
// until the true spectrum has been unfolded. Use the provided target error (error) and the maximum
// number of MLEM iterations (cutoff) to determine when to cease execution of the algorithm.
//==================================================================================================
int runMLEM(int cutoff, double error, int num_measurements, int num_bins, std::vector<double> &measurements, std::vector<double> &spectrum, std::vector<std::vector<double>> &nns_response, std::vector<double> &mlem_ratio) {
    int mlem_index; // index of MLEM iteration

    for (mlem_index = 1; mlem_index < cutoff; mlem_index++) {
        mlem_ratio.clear(); // wipe previous ratios for each iteration

        // vector that stores the MLEM-estimated data to be compared with measured data
        std::vector<double> mlem_estimate;

        // Apply system matrix, the nns_response, to current spectral estimate to get MLEM-estimated
        // data. Save results in mlem_estimate
        // Units: mlem_estimate [cps] = nns_response [cm^2] x spectru  [cps / cm^2]
        for(int i_meas = 0; i_meas < num_measurements; i_meas++)
        {
            double temp_value = 0;
            for(int i_bin = 0; i_bin < num_bins; i_bin++)
            {
                temp_value += nns_response[i_meas][i_bin]*spectrum[i_bin];
            }
            mlem_estimate.push_back(temp_value);
        }

        // Calculate ratio between each measured data point and corresponding MLEM-estimated data point
        for(int i_meas = 0; i_meas < num_measurements; i_meas++)
        {
            mlem_ratio.push_back(measurements[i_meas]/mlem_estimate[i_meas]);
        }

        // matrix that stores the correction factor to be applied to each MLEM-estimated spectral value
        std::vector<double> mlem_correction;

        // Create the correction factors to be applied to MLEM-estimated spectral values:
        //  - multiply transpose system matrix by ratio values
        for(int i_bin = 0; i_bin < num_bins; i_bin++)
        {
            double temp_value = 0;
            for(int i_meas = 0; i_meas < num_measurements; i_meas++)
            {
                temp_value += nns_response[i_meas][i_bin]*mlem_ratio[i_meas];
            }
            mlem_correction.push_back(temp_value);
        }

        // matrix that stores the normalization factors (sensitivity) to be applied to each MLEM-estimated
        // spectral value
        std::vector<double> mlem_normalization;

        // Create the normalization factors to be applied to MLEM-estimated spectral values:
        //  - each element of f stores the sum of 8 elements of the transpose system (response)
        //    matrix. The 8 elements correspond to the relative contributions of each MLEM-estimated
        //    data point to the MLEM-estimated spectral value.
        for(int i_bin = 0; i_bin < num_bins; i_bin++)
        {
            double temp_value = 0;
            for(int i_meas = 0; i_meas < num_measurements; i_meas++)
            {
                temp_value += nns_response[i_meas][i_bin];
            }
            mlem_normalization.push_back(temp_value);
        }

        // Apply correction factors and normalization to get new spectral estimate
        for(int i_bin=0; i_bin < num_bins; i_bin++)
        {
            spectrum[i_bin] = (spectrum[i_bin]*mlem_correction[i_bin]/mlem_normalization[i_bin]);
        }

        // End MLEM iterations if ratio between measured and MLEM-estimated data points is within
        // tolerace specified by 'error'
        bool continue_mlem = false;
        for (int i_meas=0; i_meas < num_measurements; i_meas++) {
            if (mlem_ratio[i_meas] >= (1+error) || mlem_ratio[i_meas] <= (1-error)) {
                continue_mlem = true;
                break;
            }   
        }
        if (!continue_mlem) {
            break;
        }
    }

    return mlem_index;
}

//==================================================================================================
// Calculate the total ambient dose equivalent rate associated with the provided spectrum using
// weighting factors (binned by energy) provided in ICRP 74.
//==================================================================================================
double calculateDose(int num_bins, std::vector<double> &spectrum, std::vector<double> &icrp_factors) {
    double ambient_dose_eq = 0;
    int s_to_hr = 3600;
    double msv_to_psv = 1e-9;

    // Sum the dose contributions from each energy bin
    for (int i_bins=0; i_bins < num_bins; i_bins++) {
        ambient_dose_eq += spectrum[i_bins]*icrp_factors[i_bins];
    }

    // Convert from [pSv/s] to [mSv/hr]
    ambient_dose_eq *= (s_to_hr)*(msv_to_psv);

    return ambient_dose_eq;
}

//==================================================================================================
// Calculate the root-mean-square deviation of sampled vectors from a "true" vector.
//==================================================================================================
int calculateRMSD_vector(int num_samples, std::vector<double> &true_vector, std::vector<std::vector<double>> &sampled_vectors, std::vector<double> &rms_differences) {
    int true_vector_size = true_vector.size();

    // Loop over each element in the true vector
    for (int i = 0; i < true_vector_size; i++) {
        double sum_sq_diff = 0;
        // Sum the square difference of poisson sampled values from the true value
        for (int i_samp = 0; i_samp < num_samples; i_samp++) {
            sum_sq_diff += ((true_vector[i] - sampled_vectors[i_samp][i])*(true_vector[i] - sampled_vectors[i_samp][i]));
        }
        double avg_sq_diff = sum_sq_diff/num_samples;
        double rms_diff = sqrt(avg_sq_diff);
        rms_differences.push_back(rms_diff);
    }
    return 1;
}

//==================================================================================================
// Calculate the root-mean-square deviation of a vector of values from a "true" value.
//==================================================================================================
double calculateRMSD(int num_samples, double true_value, std::vector<double> &sample_vector) {
    double sum_sq_diff = 0;

    // Sum the square difference of poisson sampled values from the true value
    for (int i_samp = 0; i_samp < num_samples; i_samp++) {
        sum_sq_diff += ((true_value - sample_vector[i_samp])*(true_value - sample_vector[i_samp]));
    }
    double avg_sq_diff = sum_sq_diff/num_samples;
    double rms_diff = sqrt(avg_sq_diff);

    return rms_diff;
}

//==================================================================================================
// Plot a single flux spectrum (and its uncertainty) as a function of energy. Output the generated
// plot to a file using arguments passed to the function.
//==================================================================================================
int plotSpectrum(std::string figure_file_pre, std::string figure_file_suf, std::string irradiation_conditions, int num_measurements, int num_bins, std::vector<double> &energy_bins, std::vector<double> &spectrum, std::vector<double> &spectrum_uncertainty) {
    
    // Convert vectors to arrays for input into ROOT functions
    double ini_line[num_bins];
    double bins_line[num_bins];
    double_t edges[num_bins];

    for (int i_bin = 0; i_bin < num_bins; i_bin++)
    {
        ini_line[i_bin] = spectrum[i_bin];
        bins_line[i_bin] = energy_bins[i_bin];
        edges[i_bin] = bins_line[i_bin];
    }

    // Setup plot of the spectrum
    int NBINS = num_bins-1;

    TCanvas *c1 = new TCanvas("c1","c1",800,600); // Resulution of the graph (px) specified in parameters
    TH1F *h1 = new TH1F("h1","h1",NBINS,edges);

    for (int i_bin = 0; i_bin < num_bins; i_bin++)
    {
        h1->Fill(bins_line[i_bin], ini_line[i_bin]);
    }

    h1->SetStats(0);   // Do not show the stats (mean and standard deviation);
    h1->SetLineColor(kBlue);
    h1->SetLineWidth(1);
    h1->SetTitle("NEUTRON SPECTRUM");
    h1->GetXaxis()->SetTitleOffset(1.4);
    h1->GetXaxis()->CenterTitle();
    h1->SetXTitle("Energy [MeV]");
    h1->GetYaxis()->SetTitleOffset(1.4);
    h1->GetYaxis()->CenterTitle();
    h1->SetYTitle("Fluence Rate [ncm^(-2)s^(-1)]");
    h1->Draw("HIST");  // Draw the histogram without the error bars;

    // Uncertainty plotting:

    // convert spectrum uncertainy from vector to array
    // store average bin values of adjacent bins
    double_t s_line[NBINS]; 
    double_t bins_line_avr[NBINS]; 

    for (int i_nbin = 0; i_nbin < NBINS; i_nbin++)
    {
        s_line[i_nbin] = spectrum_uncertainty[i_nbin];
        bins_line_avr[i_nbin] = (bins_line[i_nbin] + bins_line[i_nbin+1])/2;
    }

    // Setup plot of the uncertainty
    TGraphErrors *ge = new TGraphErrors(NBINS, bins_line_avr, ini_line, 0, s_line);
    ge->SetFillColor(3);
    ge->SetFillStyle(3003);
    ge->Draw("P3");

    c1->SetLogx();

    c1->Update();

    c1->Modified();

    // Output the plot to file
    std::ostringstream figure_file;
    figure_file << figure_file_pre << irradiation_conditions << figure_file_suf;
    const char *cstr_figure_file = figure_file.str().c_str();
    c1->Print(cstr_figure_file);

    return 1;
}