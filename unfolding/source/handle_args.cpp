//**************************************************************************************************
// The functions included in this module aid in processing arguments passed to the main function
// of the neutron unfolding program.
//**************************************************************************************************

#include "handle_args.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdlib.h>
#include <vector>

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
int setfile(std::vector<std::string> &arg_vector, std::string arg_string, 
    std::string default_filename, std::string &filename) 
{
    // Place an iterator at an index of the vector, where the matching arg_string was found
    std::vector<std::string>::iterator iter_args = std::find(arg_vector.begin(), arg_vector.end(), arg_string);

    // If a matching argument was found, ensure that a subsequent argument (i.e. filename to use) 
    // was provided. Apply it if so. If not, throw an error
    if( iter_args != arg_vector.end()) {
        iter_args += 1;
        if (iter_args != arg_vector.end()) {
            filename = *iter_args;
        }
        else {
            // throw error saying no file provided for *iter_measurements -= 1
            iter_args -= 1;
            throw std::logic_error("Error: no file provided for argument " + *iter_args);
        }
    }
    // If no match was found for the target arg_string within arg_vector, use the default filename
    // i.e. The user did not specify which file to use
    else {
        filename = default_filename;
    }

    return 1;
}

//==================================================================================================
// Identify any options that were provided to the main function that are not supported and notify
// the user.
// Args:
//  - arg_vector: Vector that contains all arguments passed to the main function
//  - input_file_flags: Vector that contains all valid arguments
//==================================================================================================
void checkUnknownParameters(std::vector<std::string> &arg_vector, std::vector<std::string> input_file_flags) {
    // Loop through all arguments
    for (int i=0; i<arg_vector.size(); i++) {
        // iterate through input_file_flags vector to see if any items match the current argument
        std::vector<std::string>::iterator iter_args = std::find(input_file_flags.begin(), 
            input_file_flags.end(), arg_vector[i]);
        
        // If not match was found, notify user
        if(iter_args == input_file_flags.end()) {
            throw std::logic_error("Error: Ignored unknown argument " + arg_vector[i]);
            // std::cout << "Warning: Ignored unknown argument " << arg_vector[i] << "\n";
        }
        // If a match was found, skip the next argument (should be the filename)
        else {
            i += 1;
        }
    }
}
