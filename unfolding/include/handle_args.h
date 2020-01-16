#ifndef HANDLEARGS_H
#define HANDLEARGS_H

#include <stdlib.h>
#include <fstream>
#include <vector>
#include <algorithm>

int setfile(std::vector<std::string> &arg_vector, std::string arg_string, 
    std::string default_filename, std::string &filename
);
void checkUnknownParameters(std::vector<std::string> &arg_vector, std::vector<std::string> input_file_flags);

#endif