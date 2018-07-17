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
    trend_type = "ratio";
    auto_output_path = "output/auto.csv";
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
    else if (settings_name == "trend_type")
        this->set_trend_type(settings_value);
    else if (settings_name == "auto_output_path")
        this->set_auto_output_path(settings_value);
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
void UnfoldingSettings::set_trend_type(std::string trend_type) {
    this->trend_type = trend_type;
}
void UnfoldingSettings::set_auto_output_path(std::string auto_output_path) {
    this->auto_output_path = auto_output_path;
}


//--------------------------------------------------------------------------------------------------
// Default Constructor for SpectraSettings
//--------------------------------------------------------------------------------------------------
SpectraSettings::SpectraSettings() {
    input_filename = "output_spectra.csv";
    input_dir = "output/";
    output_filename = "neutron_spectra.png";
    output_dir = "output/";
    title = "";
    x_label = "Energy (MeV)";
    y_label = "Fluence (n #upoint cm^{-2} s^{-1})";
    y_min = 0; // if plotting detects max and min are the same, use default limits
    y_max = 0;
    x_res = 800;
    y_res = 600;
    y_num_divs = 0;
    legend_entries = {};
    color_series = {"#000000","#C63822","#607FD5","#55A961"};
    color_error = {"#333333","#E79A9F","#6B8EF0","#69CF77"};
    show_error = {1};
    line_style = {1};
    line_width = {1};
    border_width = 1;
    legend_coords = {0.15,0.75,0.4,0.85};
    textbox = 0;
    textbox_coords = {0.15,0.4,0.4,0.6};
    textbox_text = {};
    plot_per_mu = 0;
    number_mu = {0};
    duration = {0};
}

// Apply a value to a setting:
void SpectraSettings::set_setting(std::string settings_name, std::string settings_value) {
    if (settings_value == ""){
        // If no setting value provided, do not apply anything
    }
    else if (settings_name == "input_filename")
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
    else if (settings_name == "legend_entries")
        this->set_legend_entries(settings_value);
    else if (settings_name == "color_series")
        this->set_color_series(settings_value);
    else if (settings_name == "color_error")
        this->set_color_error(settings_value);
    else if (settings_name == "show_error")
        this->set_show_error(settings_value);
    else if (settings_name == "line_style")
        this->set_line_style(settings_value);
    else if (settings_name == "line_width")
        this->set_line_width(settings_value);
    else if (settings_name == "border_width")
        this->set_border_width(settings_value);
    else if (settings_name == "legend_coords")
        this->set_legend_coords(settings_value);
    else if (settings_name == "textbox")
        this->set_textbox(settings_value);
    else if (settings_name == "textbox_coords")
        this->set_textbox_coords(settings_value);
    else if (settings_name == "textbox_text")
        this->set_textbox_text(settings_value);
    else if (settings_name == "plot_per_mu")
        this->set_plot_per_mu(settings_value);
    else if (settings_name == "number_mu")
        this->set_number_mu(settings_value);
    else if (settings_name == "duration")
        this->set_duration(settings_value);
    // else
    //     throw std::logic_error("Unrecognized setting: " + settings_name + ". Please refer to the README for allowed settings");
}

void SpectraSettings::set_input_filename(std::string input_filename) {
    this->input_filename = input_filename;
}
void SpectraSettings::set_input_dir(std::string input_dir) {
    this->input_dir = input_dir;
}
void SpectraSettings::set_output_filename(std::string output_filename) {
    this->output_filename = output_filename;
}
void SpectraSettings::set_output_dir(std::string output_dir) {
    this->output_dir = output_dir;
}
void SpectraSettings::set_title(std::string title) {
    this->title = title;
}
void SpectraSettings::set_x_label(std::string x_label) {
    this->x_label = x_label;
}
void SpectraSettings::set_y_label(std::string y_label) {
    this->y_label = y_label;
}
void SpectraSettings::set_y_min(std::string y_min) {
    this->y_min = stod(y_min);
}
void SpectraSettings::set_y_max(std::string y_max) {
    this->y_max = stod(y_max);
}
void SpectraSettings::set_x_res(std::string x_res) {
    this->x_res = stoi(x_res);
}
void SpectraSettings::set_y_res(std::string y_res) {
    this->y_res = stoi(y_res);
}
void SpectraSettings::set_y_num_divs(std::string y_num_divs) {
    this->y_num_divs = stoi(y_num_divs);
}
void SpectraSettings::set_legend_entries(std::string legend_entries) {
    stringToSVector(legend_entries,this->legend_entries);
}
void SpectraSettings::set_color_series(std::string color_series) {
    stringToSVector(color_series,this->color_series);
}
void SpectraSettings::set_color_error(std::string color_error) {
    stringToSVector(color_error,this->color_error);
}
void SpectraSettings::set_show_error(std::string show_error) {
    stringToIVector(show_error,this->show_error);
}
void SpectraSettings::set_line_style(std::string line_style) {
    stringToIVector(line_style,this->line_style);
}
void SpectraSettings::set_line_width(std::string line_width) {
    stringToIVector(line_width,this->line_width);
}
void SpectraSettings::set_border_width(std::string border_width) {
    this->border_width = stoi(border_width);
}
void SpectraSettings::set_legend_coords(std::string legend_coords) {
    stringToDVector(legend_coords,this->legend_coords);
}
void SpectraSettings::set_textbox(std::string textbox) {
    this->textbox = stoi(textbox);
} 
void SpectraSettings::set_textbox_coords(std::string textbox_coords) {
    stringToDVector(textbox_coords,this->textbox_coords);
}
void SpectraSettings::set_textbox_text(std::string textbox_text) {
    stringToSVector(textbox_text,this->textbox_text);
}
void SpectraSettings::set_plot_per_mu(std::string plot_per_mu) {
    this->plot_per_mu = stoi(plot_per_mu);
} 
void SpectraSettings::set_number_mu(std::string number_mu) {
    stringToIVector(number_mu,this->number_mu);
} 
void SpectraSettings::set_duration(std::string duration) {
    stringToIVector(duration,this->duration);
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
    x_log = 0;
    x_num_divs = 0;
    y_num_divs = 0;
    legend_entries = {};
    // color_series = {"#000000","#C63822","#607FD5","#55A961"};
    color_series = {"#333333","#4e79a7","#59a14f","#9c755f","#f29e2b","#edc948","#bab0ac","#e15759","#b07aa1","#76b7b2","#ff9da7"};
    // color_error = {"#333333","#E79A9F","#6B8EF0","#69CF77"};
    // show_error;
    line_style = {1};
    line_width = {2};
    border_width = 1;
    legend_coords = {0.15,0.65,0.4,0.85};
    textbox = 0;
    textbox_coords = {0.15,0.4,0.4,0.6};
    textbox_text = {};
    plot_type = {"l"};
    legend_border_size = 0;
    legend_text_size = 0.035;
    marker_style = {20};
    marker_size = {5};
    margin_left = 0.1;
    margin_right = 0.1;
    margin_top = 0.1;
    margin_bottom = 0.1;
    x_label_offset = 1;
    y_label_offset = 0;
}

// Apply a value to a setting:
void PlotSettings::set_setting(std::string settings_name, std::string settings_value) {
    if (settings_value == ""){
        // If no setting value provided, do not apply anything
    }
    else if (settings_name == "input_filename")
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
    else if (settings_name == "x_log")
        this->set_x_log(settings_value);
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
    else if (settings_name == "border_width")
        this->set_border_width(settings_value);
    else if (settings_name == "legend_coords")
        this->set_legend_coords(settings_value);
    else if (settings_name == "textbox")
        this->set_textbox(settings_value);
    else if (settings_name == "textbox_coords")
        this->set_textbox_coords(settings_value);
    else if (settings_name == "textbox_text")
        this->set_textbox_text(settings_value);
    else if (settings_name == "plot_type")
        this->set_plot_type(settings_value);
    else if (settings_name == "legend_border_size")
        this->set_legend_border_size(settings_value);
    else if (settings_name == "legend_text_size")
        this->set_legend_text_size(settings_value);
    else if (settings_name == "marker_style")
        this->set_marker_style(settings_value);
    else if (settings_name == "marker_size")
        this->set_marker_size(settings_value);
    else if (settings_name == "margin_left")
        this->set_margin_left(settings_value);
    else if (settings_name == "margin_right")
        this->set_margin_right(settings_value);
    else if (settings_name == "margin_top")
        this->set_margin_top(settings_value);
    else if (settings_name == "margin_bottom")
        this->set_margin_bottom(settings_value);
    else if (settings_name == "x_label_offset")
        this->set_x_label_offset(settings_value);
    else if (settings_name == "y_label_offset")
        this->set_y_label_offset(settings_value);
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
    this->x_min = stod(x_min);
}
void PlotSettings::set_x_max(std::string x_max) {
    this->x_max = stod(x_max);
}
void PlotSettings::set_y_min(std::string y_min) {
    this->y_min = stod(y_min);
}
void PlotSettings::set_y_max(std::string y_max) {
    this->y_max = stod(y_max);
}
void PlotSettings::set_x_res(std::string x_res) {
    this->x_res = stoi(x_res);
}
void PlotSettings::set_y_res(std::string y_res) {
    this->y_res = stoi(y_res);
}
void PlotSettings::set_x_log(std::string x_log) {
    this->x_log = stoi(x_log);
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
void PlotSettings::set_border_width(std::string border_width) {
    this->border_width = stoi(border_width);
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
void PlotSettings::set_plot_type(std::string plot_type) {
    stringToSVector(plot_type,this->plot_type);
}
void PlotSettings::set_legend_border_size(std::string legend_border_size) {
    this->legend_border_size = stoi(legend_border_size);
}
void PlotSettings::set_legend_text_size(std::string legend_text_size) {
    this->legend_text_size = stod(legend_text_size);
}
void PlotSettings::set_marker_style(std::string marker_style) {
    stringToIVector(marker_style,this->marker_style);
}
void PlotSettings::set_marker_size(std::string marker_size) {
    stringToIVector(marker_size,this->marker_size);
}
void PlotSettings::set_margin_left(std::string margin_left) {
    this->margin_left = stod(margin_left);
}
void PlotSettings::set_margin_right(std::string margin_right) {
    this->margin_right = stod(margin_right);
}
void PlotSettings::set_margin_top(std::string margin_top) {
    this->margin_top = stod(margin_top);
}
void PlotSettings::set_margin_bottom(std::string margin_bottom) {
    this->margin_bottom = stod(margin_bottom);
}
void PlotSettings::set_x_label_offset(std::string x_label_offset) {
    this->x_label_offset = stod(x_label_offset);
}
void PlotSettings::set_y_label_offset(std::string y_label_offset) {
    this->y_label_offset = stod(y_label_offset);
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
    border_width = 1;
}

// Apply a value to a setting:
void SurfaceSettings::set_setting(std::string settings_name, std::string settings_value) {
    if (settings_value == ""){
        // If no setting value provided, do not apply anything
    }
    else if (settings_name == "input_filename")
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
    else if (settings_name == "border_width")
        this->set_border_width(settings_value);

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
void SurfaceSettings::set_border_width(std::string border_width) {
    this->border_width = stoi(border_width);
}