//**************************************************************************************************
// The functions included in this module aid in reading from files and writing to them.
//**************************************************************************************************

#include "fileio.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <map>

// Constants
std::string DOSE_HEADERS[] = {
    "Irradiation Conditions",
    "Dose Rate (mSv/hr)",
    "RMS Error (mSv/hr)"
};

std::string UNCERTAINTY_SUFFIX = "_ERROR";

//==================================================================================================
// Determine if a file is empty
//==================================================================================================
bool is_empty(std::ifstream& pFile)
{
    return pFile.peek() == std::ifstream::traits_type::eof();
}


//==================================================================================================
// Retrieve settings from a configuration file 'config_file' and save values in relevant variables
// Args:
//  - config_file: filename of the settings file
//  - cutoff: assign the maximum number of MLEM iterations
//  - norm: assign the vendor-specified NNS normalization factor
//  - error: assign the maximum allowed error from unity in the ratio between measured and MLEM
//              estimated measurements
//  - f_factor: assign the conversion factor between measured current in fA and counts per second
//  - num_poisson_samples: assign the number of poisson samples to be taken for uncertainty
//              estimation
//==================================================================================================
int setSettings(std::string config_file, std::string algorithm_name, int &cutoff, double &norm, double &error, double &f_factor, double &beta, int &num_poisson_samples) {
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
    // Different algorithms require different settings
    if (algorithm_name == "mlem") {
        cutoff = (int)settings[0];
        norm = settings[1];
        error = settings[2];
        f_factor = settings[3];
        num_poisson_samples = settings[4];
    }
    else if (algorithm_name == "map") {
        cutoff = (int)settings[0];
        norm = settings[1];
        error = settings[2];
        f_factor = settings[3];
        num_poisson_samples = settings[4];
        beta = settings[5];
    }
    if (algorithm_name == "auto") {
        norm = settings[0];
        error = settings[1];
        f_factor = settings[2];
        num_poisson_samples = settings[3];
        beta = settings[4];
    }

    return true;
}

//==================================================================================================
// Retrieve settings from a configuration file 'config_file' and save values in a map
// Args:
//  - config_file: filename of the settings file
//  - settings: map variable that contains setting key:value pairs (key and value are strings)
//==================================================================================================
int setPlotSettings(std::string config_file, std::map<std::string,std::string>& settings) {
    std::ifstream cfile(config_file);
    std::string line;

    // If file is able to be read
    if (cfile.is_open())
    {
        // loop through each line in the file, extract the value for each setting into 'token'
        while ( getline (cfile,line) )
        {
            // settings format: setting_name=value
            std::string delimiter = "=";
            std::string setting_name = line.substr(0,line.find(delimiter)); // substring from 0 to '='
            std::string setting_value = line.substr(line.find(delimiter)+1); // substring from '=' to end of string
        
            bool valid_setting = checkStringMap(setting_name,settings);
            // Throw error if setting is not valid
            if (!valid_setting) {
                throw std::logic_error("Unrecognized setting: " + setting_name + ". Please refer to the README for allowed settings");
            }

            settings[setting_name] = setting_value;
        }
        cfile.close();
    }
    // Problem opening the file
    else {
        throw std::logic_error("Unable to access file: " + config_file);
    }

    return 1;
}

//==================================================================================================
// Check if a string exists within a vector of strings 
// Args:
//  - item: string to check if exists
//  - allowed_items: vector of strings to be searched
//==================================================================================================
bool checkStringVector(std::string item, std::vector<std::string>& allowed_items) {
    if(std::find(allowed_items.begin(), allowed_items.end(), item) != allowed_items.end())
    {
        return true;
    }
    else {
        return false;
    }
}

//==================================================================================================
// Check if a string key (NOTE: not a value!) exists within a map of string:string pairs
// Args:
//  - item: string to check if matches a key in the map
//  - test_map: map of strings to be searched for a particular key
//==================================================================================================
bool checkStringMap(std::string test_key, std::map<std::string, std::string>& test_map) {
    std::map<std::string, std::string>::iterator it = test_map.find(test_key);

    if(it != test_map.end())
    {
        return true;
    }
    else {
        return false;
    }
}

//==================================================================================================
// Read measurement data from file
// Args:
//  - input_file: filename of the measurement file
//  - irradiation_conditions: assign a string representing the measurement set (e.g. 10X_couch1)
//  - dose_mu: assign the delivered dose per measurement in MU
//  - doserate_mu: assign the delivered dose rate in MU/min
//  - duration: assign the measurement duration in seconds
// Returns:
//  - data_vector holds measurement values in nC. It is ordered from index:0 storing the value for 7
//      moderators to index:7 storing the value for 0 moderators
//==================================================================================================
std::vector<double> getMeasurements(std::string input_file, std::string &irradiation_conditions, double &dose_mu, double &doserate_mu, int &duration) {
    std::ifstream ifile(input_file);
    if (!ifile.is_open()) {
        //throw error
        std::cout << "Unable to open input file: " + input_file + '\n';
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
    duration = atoi(t_string.c_str());

    // Loop through file, get measurement data
    std::string line;
    std::vector<double> data_vector;
    while (getline(ifile,line)) {
        std::istringstream line_stream(line);
        std::string stoken; // store individual values between delimiters on a line

        // Loop through each line, delimiting at commas
        while (getline(line_stream, stoken, ',')) {
            // Convert negative measurements to positive
            double measurement = atof(stoken.c_str());
            if (measurement < 0) {
                measurement *= -1;
            }
            data_vector.push_back(measurement); // add data to the vector
        }
    }

    ifile.close();
    std::cout << "Data successfully retrieved from " + input_file + '\n';
    return data_vector;
}


//==================================================================================================
// Save calculated dose (and its error) to file (append to existing data in the file)
// Args:
//  - dose_file: filename to which results are saved
//  - irradiation_conditions: string that specifies measurement conditions
//  - dose: the calculated ambient dose equivalent rate
//  - dose_uncertainty: the uncertainty in the ambient dose equivalent rate
//==================================================================================================
int saveDose(std::string dose_file, std::string irradiation_conditions, double dose, double dose_uncertainty) {
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
    new_data_stream << irradiation_conditions << "," << dose << "," << dose_uncertainty << "\n";
    std::string new_data = new_data_stream.str();
    dfile << new_data; 

    dfile.close();
    return 1;
}


//==================================================================================================
// Save calculated spectrum (and uncertainty spectrum) to file (append to existing data in the file)
// Subsequent entries are appended as new rows (facilitate processing using ROOT)
//
// Args:
//  - spectrum_file: filename to which results are saved
//  - num_bins: The number of energy bins
//  - irradiation_conditions: string that specifies measurement conditions
//  - spectrum: the unfolded neutron flux spectrum
//  - spectrum_uncertainty: the uncertainty spectrum for the neutron flux
//  - energy_bins: the energy bins corresponding to spectral values
//==================================================================================================
int saveSpectrumAsRow(std::string spectrum_file, int num_bins, std::string irradiation_conditions, std::vector<double>& spectrum, std::vector<double> &spectrum_uncertainty, std::vector<double>& energy_bins) {
    // determine if file exists
    std::ifstream rfile(spectrum_file);
    bool file_empty = is_empty(rfile);
    if (!rfile.is_open()) {
        //throw error
        throw std::logic_error("Unable to access file: " + spectrum_file);
    }
    rfile.close();

    // Rewrite file using lines stored in vector 'sfile_lines'
    std::ofstream sfile(spectrum_file, std::ios::app);

    // If the file is empty, make the first line the energy bins
    if (file_empty) {
        std::ostringstream line_stream;
        line_stream << "Energy (MeV)" << ",";
        for (int i_bin = 0; i_bin < num_bins; i_bin++) {
            line_stream << energy_bins[i_bin];
            if (i_bin < num_bins-1) {
                line_stream << ",";
            }
        }
        sfile << line_stream.str();
    }
    
    // Append lines for the unfolded spectrum & its uncertainty
    std::ostringstream spectrum_stream;
    std::ostringstream error_stream;

    spectrum_stream << irradiation_conditions << ",";
    error_stream << irradiation_conditions << UNCERTAINTY_SUFFIX << ",";

    for (int i_bin = 0; i_bin < num_bins; i_bin++) {
        spectrum_stream << spectrum[i_bin];
        error_stream << spectrum_uncertainty[i_bin];
        if (i_bin < num_bins-1) {
            spectrum_stream << ",";
            error_stream << ",";
        }
    }

    sfile << "\r\n" << spectrum_stream.str();
    sfile << "\r\n" << error_stream.str();    

    return 1;
}

//==================================================================================================
// Save calculated spectrum (and uncertainty spectrum) to file (append to existing data in the file)
// Subsequent entries are saved as new columns (facilitate import into spreadsheet)
//
// Args:
//  - spectrum_file: filename to which results are saved
//  - irradiation_conditions: string that specifies measurement conditions
//  - spectrum: the unfolded neutron flux spectrum
//  - spectrum_uncertainty: the uncertainty spectrum for the neutron flux
//  - energy_bins: the energy bins corresponding to spectral values
//==================================================================================================
int saveSpectrumAsColumn(std::string spectrum_file, std::string irradiation_conditions, std::vector<double>& spectrum, std::vector<double> &spectrum_uncertainty, std::vector<double>& energy_bins) {
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
// Read a 1D CSV file (i.e. one value per row) and store the data in a vector
// Args:
//  - file_name: the name of the file to be read
//  - input_vector: the vector that will be assigned values read in from the input file (note passed
//      by reference)
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
// Args:
//  - file_name: the name of the file to be read
//  - input_vector: the vector that will be assigned values read in from the input file (note passed
//      by reference)
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
// Read a CSV file containing neutron spectra, formatted as follows:
//  - The first row contains the energy bins
//  - Subsequent rows contain alternating spectra and their uncertainty spectra
//  - The first column of every row contains a header for that row (non-numeric)
//
// Args:
//  - file_name: the name of the file to be read
//  - header_vector: the vector that will be assigned header values
//  - energy_bins: the vector that will house the energy bins
//  - spectra_vector: the vector that will be assigned spectra values from the input file
//  - error_vector: the vector that will be assigned uncertainties from the input file
//==================================================================================================
int readSpectra(std::string file_name, std::vector<std::string>& header_vector, std::vector<double>& energy_bins, std::vector<std::vector<double>>& spectra_vector, std::vector<std::vector<double>>& error_vector) {
    std::ifstream ifile(file_name);
    std::string iline;

    if (!ifile.is_open()) {
        //throw error
        throw std::logic_error("Unable to open input file: " + file_name);
    }

    // Loop through each line in the file
    int i_row = 0;
    while (getline(ifile,iline)) {
        std::vector<double> new_row;
        std::istringstream line_stream(iline);
        std::string stoken; // store individual values between delimiters on a line

        // Loop through the line, delimiting at commas
        int i_col = 0;
        while (getline(line_stream, stoken, ',')) {
            // The first token is the header
            if (i_col == 0) {
                if (i_row % 2 == 1) {
                    header_vector.push_back(stoken);
                }
            }
            else {
                new_row.push_back(atof(stoken.c_str())); // add data to the vector
            }
            i_col += 1;
        }

        // The first row contains the energy bins
        if (i_row == 0) {
            energy_bins = new_row;
        }
        // Remaining rows alternate between a spectrum and an uncertainty spectrum
        else if (i_row % 2 == 1) {
            spectra_vector.push_back(new_row);
        }
        else {
            error_vector.push_back(new_row);
        }
        i_row += 1;
    }

    return 1;
}

//==================================================================================================
// Compare the size of two vectors (indicated by reference_size and test_size). If the size of test
// does not match the size of reference, throw an error with a message constructed using text passed
// as function parameters (reference_string and test_string)
// Args:
//  - reference_size: length of the reference vector
//  - reference_string: string to be included in the error message that represents the reference
//      vector
//  - test_size: length of the test vector
//  - test_string: string to be included in the error message that represents the test
//      vector
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
// Generate a textfile report of pertinent information from execution of this program. The contents
// of this function are separated by headers indicating the type of information printed to the
// report in the the corresponding section.
//==================================================================================================
int prepareReport(std::string report_file, std::string irradiation_conditions, std::vector<std::string> &input_files, std::vector<std::string> &input_file_flags, std::string algorithm_name, int cutoff, double error, double norm, double f_factor, double beta, int num_measurements, int num_bins, int num_poisson_samples, std::vector<double>& measurements_nc, double dose_mu, double doserate_mu, int duration, std::vector<double>& energy_bins, std::vector<double>& initial_spectrum, std::vector<std::vector<double>>& nns_response, int num_iterations, std::vector<double>& mlem_ratio, double dose, double s_dose, std::vector<double>& spectrum, std::vector<double>& spectrum_uncertainty, std::vector<double>& icrp_factors, std::string git_commit) {
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
    rfile << std::left << std::setw(sw) << "Git commit number: " << git_commit << "\n";
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
        rfile << std::left << std::setw(cw) << num_measurements-1-i << measurements_nc[i] << "\n";
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
    rfile << "Unfolding information\n\n";
    rfile << std::left << std::setw(sw) << "Algorithm: " << algorithm_name << "\n";
    if (algorithm_name == "map") {
        rfile << std::left << std::setw(sw) << "MAP beta value: " << beta << "\n";
    }
    rfile << std::left << std::setw(sw) << "# of iterations: " << num_iterations << "/" << cutoff << "\n\n";
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
