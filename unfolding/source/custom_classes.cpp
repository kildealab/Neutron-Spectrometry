#include "custom_classes.h"
#include "fileio.h"

#include <stdlib.h>
#include <string>
#include <iostream>

//--------------------------------------------------------------------------------------------------
// Default constructor for UnfoldingSettings
//--------------------------------------------------------------------------------------------------
UnfoldingSettings::UnfoldingSettings() {
    norm = 1.21;
    error = 0.0;
    f_factor = 7.0;
    cutoff = 10000;
    num_poisson_samples = 50;
    // MAP specific
    beta = 0.0;
    prior = "mrp";
    // Optimize specific
    min_num_iterations = 500;
    max_num_iterations = 10000;
    iteration_increment = 50;
    min_beta = 0.0;
    max_beta = 1.0;
    parameter_of_interest = "fluence";
    algorithm = "mlem";
}

// Apply a value to a setting:
void UnfoldingSettings::set_setting(std::string settings_name, std::string settings_value) {
    if (settings_name == "nns_normalization")
        this->set_norm(atof(settings_value.c_str()));
    else if (settings_name == "mlem_max_error")
        this->set_error(atof(settings_value.c_str()));
    else if (settings_name == "f_factor")
        this->set_f_factor(atof(settings_value.c_str()));
    else if (settings_name == "mlem_cutoff")
        this->set_cutoff(atoi(settings_value.c_str()));
    else if (settings_name == "num_poisson_samples")
        this->set_num_poisson_samples(atoi(settings_value.c_str()));
    else if (settings_name == "beta")
        this->set_beta(atof(settings_value.c_str()));
    else if (settings_name == "prior")
        this->set_prior(settings_value);
    else if (settings_name == "min_num_iterations")
        this->set_min_num_iterations(atoi(settings_value.c_str()));
    else if (settings_name == "max_num_iterations")
        this->set_max_num_iterations(atoi(settings_value.c_str()));
    else if (settings_name == "iteration_increment")
        this->set_iteration_increment(atoi(settings_value.c_str()));
    else if (settings_name == "min_beta")
        this->set_min_beta(atof(settings_value.c_str()));
    else if (settings_name == "max_beta")
        this->set_max_beta(atof(settings_value.c_str()));
    else if (settings_name == "parameter_of_interest")
        this->set_parameter_of_interest(settings_value);
    else if (settings_name == "algorithm")
        this->set_algorithm(settings_value);
    else
        throw std::logic_error("Unrecognized setting: " + settings_name + ". Please refer to the README for allowed settings");

}

// Setter functions
void UnfoldingSettings::set_norm(double norm) {
    this->norm = norm;
}
void UnfoldingSettings::set_error(double error) {
    this->error = error;
}
void UnfoldingSettings::set_f_factor(double f_factor) {
    this->f_factor = f_factor;
}
void UnfoldingSettings::set_cutoff(int cutoff) {
    this->cutoff = cutoff;
}
void UnfoldingSettings::set_num_poisson_samples(int num_poisson_samples) {
    this->num_poisson_samples = num_poisson_samples;
}
void UnfoldingSettings::set_beta(double beta) {
    this->beta = beta;
}
void UnfoldingSettings::set_prior(std::string prior) {
    this->prior = prior;
}
void UnfoldingSettings::set_min_num_iterations(int min_num_iterations) {
    this->min_num_iterations = min_num_iterations;
}
void UnfoldingSettings::set_max_num_iterations(int max_num_iterations) {
    this->max_num_iterations = max_num_iterations;
}
void UnfoldingSettings::set_iteration_increment(int iteration_increment) {
    this->iteration_increment = iteration_increment;
}
void UnfoldingSettings::set_min_beta(double min_beta) {
    this->min_beta = min_beta;
}
void UnfoldingSettings::set_max_beta(double max_beta) {
    this->max_beta = max_beta;
}
void UnfoldingSettings::set_parameter_of_interest(std::string parameter_of_interest) {
    this->parameter_of_interest = parameter_of_interest;
}
void UnfoldingSettings::set_algorithm(std::string algorithm) {
    this->algorithm = algorithm;
}

//--------------------------------------------------------------------------------------------------
// Default Constructor for PlotSettings
//--------------------------------------------------------------------------------------------------
PlotSettings::PlotSettings() {
    input_filename = "poi_output_mlem.csv";
    input_dir = "output/";
    output_filename = "poi_output_mlem.png";
    output_dir = "output/";
    title = "";
    x_label = "";
    y_label = "";
    y_min = 0; // if plotting detects max and min are the same, use default limits
    y_max = 0;
    x_min = 0; // if plotting detects max and min are the same, use default limits
    x_max = 0;
    x_res = 800;
    y_res = 600;
    x_num_divs = 0;
    y_num_divs = 0;
    // legend_entries;
    color_series = {"#000000","#C63822","#607FD5","#55A961"};
    // color_error = {"#333333","#E79A9F","#6B8EF0","#69CF77"};
    // show_error;
    line_style = {1,1,1,1};
    line_width = {2,2,2,2};
    legend_coords = {0.15,0.65,0.4,0.85};
    textbox = 0;
    textbox_coords = {0.15,0.4,0.4,0.6};
    // textbox_text;
}

// Apply a value to a setting:
void PlotSettings::set_setting(std::string settings_name, std::string settings_value) {
    if (settings_name == "input_filename")
        this->set_input_filename(settings_value);
    else if (settings_name == "input_dir")
        this->set_input_dir(settings_value);
    else if (settings_name == "output_filename")
        this->set_output_filename(settings_value);
    else if (settings_name == "output_dir")
        this->set_output_dir(settings_value);
    else if (settings_name == "title")
        this->set_title(settings_value);
    else if (settings_name == "x_label")
        this->set_x_label(settings_value);
    else if (settings_name == "y_label")
        this->set_y_label(settings_value);
    else if (settings_name == "x_min")
        this->set_x_min(settings_value);
    else if (settings_name == "x_max")
        this->set_x_max(settings_value);
    else if (settings_name == "y_min")
        this->set_y_min(settings_value);
    else if (settings_name == "y_max")
        this->set_y_max(settings_value);
    else if (settings_name == "x_res")
        this->set_x_res(settings_value);
    else if (settings_name == "y_res")
        this->set_y_res(settings_value);
    else if (settings_name == "y_num_divs")
        this->set_y_num_divs(settings_value);
    else if (settings_name == "x_num_divs")
        this->set_x_num_divs(settings_value);
    else if (settings_name == "legend_entries")
        this->set_legend_entries(settings_value);
    else if (settings_name == "color_series")
        this->set_color_series(settings_value);
    // else if (settings_name == "color_error")
    //     this->set_color_error(settings_value);
    // else if (settings_name == "show_error")
    //     this->set_show_error(settings_value);
    else if (settings_name == "line_style")
        this->set_line_style(settings_value);
    else if (settings_name == "line_width")
        this->set_line_width(settings_value);
    else if (settings_name == "legend_coords")
        this->set_legend_coords(settings_value);
    else if (settings_name == "textbox")
        this->set_textbox(settings_value);
    else if (settings_name == "textbox_coords")
        this->set_textbox_coords(settings_value);
    else if (settings_name == "textbox_text")
        this->set_textbox_text(settings_value);
    else
        throw std::logic_error("Unrecognized setting: " + settings_name + ". Please refer to the README for allowed settings");
}

// Setter functions
void PlotSettings::set_input_filename(std::string input_filename) {
    this->input_filename = input_filename;
}
void PlotSettings::set_input_dir(std::string input_dir) {
    this->input_dir = input_dir;
}
void PlotSettings::set_output_filename(std::string output_filename) {
    this->output_filename = output_filename;
}
void PlotSettings::set_output_dir(std::string output_dir) {
    this->output_dir = output_dir;
}
void PlotSettings::set_title(std::string title) {
    this->title = title;
}
void PlotSettings::set_x_label(std::string x_label) {
    this->x_label = x_label;
}
void PlotSettings::set_y_label(std::string y_label) {
    this->y_label = y_label;
}
void PlotSettings::set_x_min(std::string x_min) {
    this->x_min = stoi(x_min);
}
void PlotSettings::set_x_max(std::string x_max) {
    this->x_max = stoi(x_max);
}
void PlotSettings::set_y_min(std::string y_min) {
    this->y_min = stoi(y_min);
}
void PlotSettings::set_y_max(std::string y_max) {
    this->y_max = stoi(y_max);
}
void PlotSettings::set_x_res(std::string x_res) {
    this->x_res = stoi(x_res);
}
void PlotSettings::set_y_res(std::string y_res) {
    this->y_res = stoi(y_res);
}
void PlotSettings::set_x_num_divs(std::string x_num_divs) {
    this->x_num_divs = stoi(x_num_divs);
}
void PlotSettings::set_y_num_divs(std::string y_num_divs) {
    this->y_num_divs = stoi(y_num_divs);
}
void PlotSettings::set_legend_entries(std::string legend_entries) {
    stringToSVector(legend_entries,this->legend_entries);
}
void PlotSettings::set_color_series(std::string color_series) {
    stringToSVector(color_series,this->color_series);
}
// void PlotSettings::set_color_error(std::string color_error) {
//     stringToSVector(color_error,this->color_error);
// }
// void PlotSettings::set_show_error(std::string show_error) {
//     stringToIVector(show_error,this->show_error);
// }
void PlotSettings::set_line_style(std::string line_style) {
    stringToIVector(line_style,this->line_style);
}
void PlotSettings::set_line_width(std::string line_width) {
    stringToIVector(line_width,this->line_width);
}
void PlotSettings::set_legend_coords(std::string legend_coords) {
    stringToDVector(legend_coords,this->legend_coords);
}
void PlotSettings::set_textbox(std::string textbox) {
    this->textbox = stoi(textbox);
} 
void PlotSettings::set_textbox_coords(std::string textbox_coords) {
    stringToDVector(textbox_coords,this->textbox_coords);
}
void PlotSettings::set_textbox_text(std::string textbox_text) {
    stringToSVector(textbox_text,this->textbox_text);
}



//--------------------------------------------------------------------------------------------------
// Default Constructor for SurfaceSettings
//--------------------------------------------------------------------------------------------------
SurfaceSettings::SurfaceSettings() {
    input_filename = "poi_output_map.csv";
    input_dir = "output/";
    output_filename = "poi_output_map.png";
    output_dir = "output/";
    title = "";
    x_label = "";
    y_label = "";
    z_label = "";
    y_min = 0; // if plotting detects max and min are the same, use default limits
    y_max = 0;
    x_min = 0; // if plotting detects max and min are the same, use default limits
    x_max = 0;
    z_min = 0; // if plotting detects max and min are the same, use default limits
    z_max = 0;
    x_res = 800;
    y_res = 600;
    x_num_divs = 206;
    z_num_divs = 203;
    color_palette = 55;
    num_color_bins = 40;
}

// Apply a value to a setting:
void SurfaceSettings::set_setting(std::string settings_name, std::string settings_value) {
    if (settings_name == "input_filename")
        this->set_input_filename(settings_value);
    else if (settings_name == "input_dir")
        this->set_input_dir(settings_value);
    else if (settings_name == "output_filename")
        this->set_output_filename(settings_value);
    else if (settings_name == "output_dir")
        this->set_output_dir(settings_value);
    else if (settings_name == "title")
        this->set_title(settings_value);
    else if (settings_name == "x_label")
        this->set_x_label(settings_value);
    else if (settings_name == "y_label")
        this->set_y_label(settings_value);
    else if (settings_name == "z_label")
        this->set_z_label(settings_value);
    else if (settings_name == "x_min")
        this->set_x_min(settings_value);
    else if (settings_name == "x_max")
        this->set_x_max(settings_value);
    else if (settings_name == "y_min")
        this->set_y_min(settings_value);
    else if (settings_name == "y_max")
        this->set_y_max(settings_value);
    else if (settings_name == "z_min")
        this->set_z_min(settings_value);
    else if (settings_name == "z_max")
        this->set_z_max(settings_value);
    else if (settings_name == "x_res")
        this->set_x_res(settings_value);
    else if (settings_name == "y_res")
        this->set_y_res(settings_value);
    else if (settings_name == "z_num_divs")
        this->set_z_num_divs(settings_value);
    else if (settings_name == "x_num_divs")
        this->set_x_num_divs(settings_value);
    else if (settings_name == "color_palette")
        this->set_color_palette(settings_value);
    else if (settings_name == "num_color_bins")
        this->set_num_color_bins(settings_value);

    else
        throw std::logic_error("Unrecognized setting: " + settings_name + ". Please refer to the README for allowed settings");
}

// Setter functions
void SurfaceSettings::set_input_filename(std::string input_filename) {
    this->input_filename = input_filename;
}
void SurfaceSettings::set_input_dir(std::string input_dir) {
    this->input_dir = input_dir;
}
void SurfaceSettings::set_output_filename(std::string output_filename) {
    this->output_filename = output_filename;
}
void SurfaceSettings::set_output_dir(std::string output_dir) {
    this->output_dir = output_dir;
}
void SurfaceSettings::set_title(std::string title) {
    this->title = title;
}
void SurfaceSettings::set_x_label(std::string x_label) {
    this->x_label = x_label;
}
void SurfaceSettings::set_y_label(std::string y_label) {
    this->y_label = y_label;
}
void SurfaceSettings::set_z_label(std::string z_label) {
    this->z_label = z_label;
}
void SurfaceSettings::set_x_min(std::string x_min) {
    this->x_min = stoi(x_min);
}
void SurfaceSettings::set_x_max(std::string x_max) {
    this->x_max = stoi(x_max);
}
void SurfaceSettings::set_y_min(std::string y_min) {
    this->y_min = stod(y_min);
}
void SurfaceSettings::set_y_max(std::string y_max) {
    this->y_max = stod(y_max);
}
void SurfaceSettings::set_z_min(std::string z_min) {
    this->z_min = stod(z_min);
}
void SurfaceSettings::set_z_max(std::string z_max) {
    this->z_max = stod(z_max);
}
void SurfaceSettings::set_x_res(std::string x_res) {
    this->x_res = stoi(x_res);
}
void SurfaceSettings::set_y_res(std::string y_res) {
    this->y_res = stoi(y_res);
}
void SurfaceSettings::set_x_num_divs(std::string x_num_divs) {
    this->x_num_divs = stoi(x_num_divs);
}
void SurfaceSettings::set_z_num_divs(std::string z_num_divs) {
    this->z_num_divs = stoi(z_num_divs);
}
void SurfaceSettings::set_color_palette(std::string color_palette) {
    this->color_palette = stoi(color_palette);
}
void SurfaceSettings::set_num_color_bins(std::string num_color_bins) {
    this->num_color_bins = stoi(num_color_bins);
}