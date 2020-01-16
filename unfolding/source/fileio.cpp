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
//  - settings: UnfoldingSettings object that will store necessary parameters input from config_file
//==================================================================================================
int setSettings(std::string config_file, UnfoldingSettings &settings) {
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
            std::string settings_name = line.substr(0,line.find(delimiter)); // substring from 0 to '='
            std::string settings_value = line.substr(line.find(delimiter)+1); // substring from '=' to end of string
            
            // Apply the setting value to the setting name
            settings.set_setting(settings_name,settings_value);
        }
        cfile.close();
    }

    // Problem opening the file
    else {
        throw std::logic_error("Unable to access configuration file: " + config_file);
    }

    return 1;
}


//==================================================================================================
// Retrieve settings from a configuration file 'config_file' and save values in a settings object
// Args:
//  - config_file: filename of the settings file
//  - settings: SpectraSettings object that will store necessary parameters input from config_file
//==================================================================================================
int setSpectraSettings(std::string config_file, SpectraSettings &settings) {
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
            std::string settings_name = line.substr(0,line.find(delimiter)); // substring from 0 to '='
            std::string settings_value = line.substr(line.find(delimiter)+1); // substring from '=' to end of string

            // Apply the setting value to the setting name
            settings.set_setting(settings_name,settings_value);
        }
        cfile.close();
    }

    // Problem opening the file
    else {
        throw std::logic_error("Unable to access configuration file: " + config_file);
    }

    return 1;
}


//==================================================================================================
// Retrieve settings from a configuration file 'config_file' and save values in a settings object
// Args:
//  - config_file: filename of the settings file
//  - settings: PlotSettings object that will store necessary parameters input from config_file
//==================================================================================================
int setPlotSettings(std::string config_file, PlotSettings &settings) {
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
            std::string settings_name = line.substr(0,line.find(delimiter)); // substring from 0 to '='
            std::string settings_value = line.substr(line.find(delimiter)+1); // substring from '=' to end of string

            // Apply the setting value to the setting name
            settings.set_setting(settings_name,settings_value);
        }
        cfile.close();
    }

    // Problem opening the file
    else {
        throw std::logic_error("Unable to access configuration file: " + config_file);
    }

    return 1;
}

//==================================================================================================
// Retrieve settings from a configuration file 'config_file' and save values in a settings object
// Args:
//  - config_file: filename of the settings file
//  - settings: SurfaceSettings object that will store necessary parameters input from config_file
//==================================================================================================
int setSurfaceSettings(std::string config_file, SurfaceSettings &settings) {
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
            std::string settings_name = line.substr(0,line.find(delimiter)); // substring from 0 to '='
            std::string settings_value = line.substr(line.find(delimiter)+1); // substring from '=' to end of string

            // Apply the setting value to the setting name
            settings.set_setting(settings_name,settings_value);
        }
        cfile.close();
    }

    // Problem opening the file
    else {
        throw std::logic_error("Unable to access configuration file: " + config_file);
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
std::vector<double> getMeasurements(UnfoldingSettings &settings) 
{
    std::ifstream ifile(settings.path_measurements);
    if (!ifile.is_open()) {
        //throw error
        throw std::logic_error("Unable to access open measurement file: " + settings.path_measurements);
    }

    // Load header information from 'ifile'
    std::string temp_irradiation_conditions;
    getline(ifile,temp_irradiation_conditions);
    // removes: carriage return '\r' from the string (which causes weird string overwriting)
    temp_irradiation_conditions.erase( std::remove(temp_irradiation_conditions.begin(), 
        temp_irradiation_conditions.end(), '\r'), temp_irradiation_conditions.end() );
    settings.irradiation_conditions = temp_irradiation_conditions;

    // If measurements in nC, extract dose, doserate, and duration information
    if (settings.meas_units == "nc") {
        std::string dose_string;
        getline(ifile,dose_string);
        settings.dose_mu = atoi(dose_string.c_str());

        std::string doserate_string;
        getline(ifile,doserate_string);
        settings.doserate_mu = atoi(doserate_string.c_str());

        std::string t_string;
        getline(ifile,t_string);
        settings.duration = atoi(t_string.c_str());
    }

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
    std::cout << "Measurements successfully retrieved from " + settings.path_measurements + '\n';
    return data_vector;
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
//  - error_lower: the lower uncertainty spectrum for the neutron flux
//  - error_upper: the upper uncertainty spectrum for the neutron flux
//  - energy_bins: the energy bins corresponding to spectral values
//==================================================================================================
int saveSpectrumAsRow(std::string spectrum_file, int num_bins, std::string irradiation_conditions, 
    std::vector<double>& spectrum, std::vector<double> &error_lower, std::vector<double> &error_upper,
    std::vector<double>& energy_bins) 
{
    std::string error_lower_label = "Lower uncertainty";
    std::string error_upper_label = "Upper uncertainty";

    // determine if file exists
    std::ifstream rfile(spectrum_file);
    bool file_empty = is_empty(rfile);
    // bool file_exists = rfile.good();
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
    std::ostringstream error_lower_stream;
    std::ostringstream error_upper_stream;

    spectrum_stream << irradiation_conditions << ",";
    error_lower_stream << error_lower_label << ",";
    error_upper_stream << error_upper_label << ",";

    for (int i_bin = 0; i_bin < num_bins; i_bin++) {
        spectrum_stream << spectrum[i_bin];
        error_lower_stream << error_lower[i_bin];
        error_upper_stream << error_upper[i_bin];
        if (i_bin < num_bins-1) {
            spectrum_stream << ",";
            error_lower_stream << ",";
            error_upper_stream << ",";
        }
    }

    sfile << "\r\n" << spectrum_stream.str();
    sfile << "\r\n" << error_lower_stream.str();   
    sfile << "\r\n" << error_upper_stream.str();    

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
int saveSpectrumAsColumn(std::string spectrum_file, std::string irradiation_conditions, std::vector<double>& spectrum, 
    std::vector<double> &spectrum_uncertainty, std::vector<double>& energy_bins) 
{
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
int readSpectra(std::string file_name, std::vector<std::string>& header_vector, std::vector<double>& energy_bins, 
    std::vector<std::vector<double>>& spectra_vector, std::vector<std::vector<double>>& error_lower_vector, 
    std::vector<std::vector<double>>& error_upper_vector, bool plot_per_mu, std::vector<int>& number_mu, 
    std::vector<int>& duration, int rows_per_spectrum) 
{
    std::ifstream ifile(file_name);
    std::string iline;

    if (!ifile.is_open()) {
        //throw error
        throw std::logic_error("Unable to open spectrum file: " + file_name);
    }

    // Loop through each line in the file
    int i_row = 0;
    int i_group = 0; // index that is incremented after processing each spectrum & its uncertainties
    while (getline(ifile,iline)) {
        std::vector<double> new_row;
        std::istringstream line_stream(iline);
        std::string stoken; // store individual values between delimiters on a line

        // Loop through the line, delimiting at commas
        int i_col = 0;
        while (getline(line_stream, stoken, ',')) {
            // The first token is the header
            if (i_col == 0) {
                if (i_row % rows_per_spectrum == 1) {
                    header_vector.push_back(stoken);
                }
            }
            // Add data to the vector
            else {
                // If plotting per MU
                if (plot_per_mu){
                    //1st row is the energies
                    if (i_row !=0) {
                        new_row.push_back(
                            atof(stoken.c_str())*duration[i_group%duration.size()]/number_mu[i_group%number_mu.size()]
                        );
                    }
                    else {
                        new_row.push_back(atof(stoken.c_str()));
                    }
                }
                // If plotting per second
                else {
                    new_row.push_back(atof(stoken.c_str()));
                }
            }
            i_col += 1;
        }

        //
        if (rows_per_spectrum == 3) {
            // The first row contains the energy bins
            if (i_row == 0) {
                energy_bins = new_row;
            }
            // Remaining rows alternate between a spectrum and upper & lower uncertainties
            else if (i_row % rows_per_spectrum == 1) {
                spectra_vector.push_back(new_row);
            }
            else if (i_row % rows_per_spectrum == 2) {
                error_upper_vector.push_back(new_row);
                i_group += 1;
            }
            else {
                error_lower_vector.push_back(new_row);
            }
        }
        else if (rows_per_spectrum == 2) {
            // The first row contains the energy bins
            if (i_row == 0) {
                energy_bins = new_row;
            }
            // Remaining rows alternate between a spectrum and the uncertainty
            else if (i_row % rows_per_spectrum == 1) {
                spectra_vector.push_back(new_row);
            }
            else {
                error_upper_vector.push_back(new_row);
                error_lower_vector.push_back(new_row);
                i_group += 1;
            }
        }
        else {
            throw std::logic_error("Incompatible number for rows_per_spectrum: " + std::to_string(rows_per_spectrum));
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
        error_message << "File dimensions mismatch: " << test_string << " (" << test_size << ") does not match " 
            << reference_string << " (" << reference_size << ")";
        throw std::logic_error(error_message.str());   
    }
    return 1;
}


//==================================================================================================
// Convert a comma-delimited string into a vector of strings.
//
// Args:
//  - test_string: the comma-delimited string to be processed
//  - result_vector: the vector that will be assigned string values
//==================================================================================================
void stringToSVector(std::string test_string, std::vector<std::string>& result_vector) {
    result_vector.clear();
    std::istringstream line_stream(test_string);
    std::string stoken; // store individual values between delimiters on a line

    // Loop through each line, delimiting at commas
    while (getline(line_stream, stoken, ',')) {
        result_vector.push_back(stoken);
    }
}

//==================================================================================================
// Convert a comma-delimited string into a vector of integers.
//
// Args:
//  - test_string: the comma-delimited string to be processed
//  - result_vector: the vector that will be assigned integer values
//==================================================================================================
void stringToIVector(std::string test_string, std::vector<int>& result_vector) {
    result_vector.clear();
    std::istringstream line_stream(test_string);
    std::string stoken; // store individual values between delimiters on a line

    // Loop through each line, delimiting at commas
    while (getline(line_stream, stoken, ',')) {
        result_vector.push_back(atoi(stoken.c_str()));
    }
}

//==================================================================================================
// Convert a comma-delimited string into a vector of floats.
//
// Args:
//  - test_string: the comma-delimited string to be processed
//  - result_vector: the vector that will be assigned integer values
//==================================================================================================
void stringToDVector(std::string test_string, std::vector<float>& result_vector) {
    result_vector.clear();
    std::istringstream line_stream(test_string);
    std::string stoken; // store individual values between delimiters on a line

    // Loop through each line, delimiting at commas
    while (getline(line_stream, stoken, ',')) {
        result_vector.push_back((double) atof(stoken.c_str()));
    }
}

//==================================================================================================
// Read a CSV file containing XYY data to be plotted that is formatted as follows
//  - The first row contains the x data
//  - Subsequent rows contain y data
//  - The first column of every row contains a header for that row (non-numeric)
//
// Args:
//  - file_name: the name of the file to be read
//  - header_vector: the vector that will be assigned header values
//  - x_data: the 1D vector that will store the x data (first row)
//  - y_data: the 2D vector that will store the y data (subsequent rows)
//==================================================================================================
int readXYYCSV(std::string file_name, std::vector<std::string>& header_vector,
    std::vector<std::vector<double>>& x_data, std::vector<std::vector<double>>& y_data) 
{
    std::ifstream ifile(file_name);
    std::string iline;

    if (!ifile.is_open()) {
        //throw error
        throw std::logic_error("Unable to open data file: " + file_name);
    }

    std::vector<double> x_row;

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
                // Only add y-row headers to the header vector
                if (i_row > 0) {
                    header_vector.push_back(stoken);
                }
            }
            // Add data to the vector
            else {
                new_row.push_back(atof(stoken.c_str()));
            }
            i_col += 1;
        }

        // The first row contains the x data
        if (i_row == 0) {
            x_row = new_row;
            x_data.push_back(x_row);
        }
        // The second row contain the y data
        else if (i_row == 1) {
            y_data.push_back(new_row);
        }
        // All remaining rows contain y data
        else {
            x_data.push_back(x_row);
            y_data.push_back(new_row);
        }
        i_row += 1;
    }
    return 1;
}

//==================================================================================================
// Read a CSV file containing XYXY data to be plotted that is formatted as follows
//  - Even rows (starting at 0) contain x data
//  - Odd rows (starting at 1) contain y data
//  - The first column of every row contains a header for that row (non-numeric)
//
// Args:
//  - file_name: the name of the file to be read
//  - header_vector: the vector that will be assigned header values
//  - x_data: the 2D vector that will store the x data (first row)
//  - y_data: the 2D vector that will store the y data (subsequent rows)
//==================================================================================================
int readXYXYCSV(std::string file_name, std::vector<std::string>& header_vector, 
    std::vector<std::vector<double>>& x_data, std::vector<std::vector<double>>& y_data) 
{
    std::ifstream ifile(file_name);
    std::string iline;

    if (!ifile.is_open()) {
        //throw error
        throw std::logic_error("Unable to open data file: " + file_name);
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
                // Only add y-row headers to the header vector
                if (i_row % 2 == 1) {
                    header_vector.push_back(stoken);
                }
            }
            // Add data to the vector
            else {
                new_row.push_back(atof(stoken.c_str()));
            }
            i_col += 1;
        }

        // Even rows (starting at 0) contain the x data
        if (i_row % 2 == 0) {
            x_data.push_back(new_row);
        }
        // Odd rows (starting at 1) contain the y data
        else {
            y_data.push_back(new_row);
        }
        i_row += 1;
    }
    return 1;
}