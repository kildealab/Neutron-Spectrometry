#ifndef CUSTOM_CLASSES_H
#define CUSTOM_CLASSES_H

#include <stdlib.h>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>

#include "physics_calculations.h"

// class Settings {
//     public:
//         Settings(){};
//         virtual ~Settings() {};
//         virtual void set_setting (std::string,std::string) = 0;
// };

class UnfoldingSettings{
// class UnfoldingSettings : public Settings{
    public:
        double norm; 
        double error; 
        double f_factor; 
        int cutoff; 
        std::string uncertainty_type;
        int num_uncertainty_samples;
        int num_meas_per_shell; 
        std::string meas_units; 

        int dose_mu;
        int doserate_mu;
        int duration;
        std::string irradiation_conditions;

        std::string path_measurements;
        std::string path_input_spectrum;
        std::string path_energy_bins;
        std::string path_system_response;
        std::string path_icrp_factors;

        std::string path_output_spectra;
        int generate_report;
        std::string path_report;
        int generate_figure;
        std::string path_figure;

        // MAP specific
        double beta; 
        std::string prior;
        //MLEM-STOP specific
        int cps_crossover;
        double sigma_j;
        // Optimize specific
        int iteration_min;
        int iteration_max;
        int iteration_increment;
        double beta_min;
        double beta_max;
        std::string parameter_of_interest;
        std::string algorithm;
        std::string trend_type;
        int derivatives;
        std::string path_output_trend;
        std::string path_ref_spectrum;

        UnfoldingSettings(); 

        void set_setting(std::string,std::string);

        void set_norm(double);
        void set_error(double);
        void set_f_factor(double);
        void set_cutoff(int);
        void set_uncertainty_type(std::string);
        void set_num_uncertainty_samples(int);
        void set_num_meas_per_shell(int);
        void set_meas_units(std::string);
        void set_dose_mu(int);
        void set_doserate_mu(int);
        void set_duration(int);
        void set_irradiation_conditions(std::string);
        void set_beta(double);
        void set_prior(std::string);
        void set_cps_crossover(int);
        void set_sigma_j(double);
        void set_iteration_min(int);
        void set_iteration_max(int);
        void set_iteration_increment(int);
        void set_beta_min(double);
        void set_beta_max(double);
        void set_parameter_of_interest(std::string);
        void set_algorithm(std::string);
        void set_trend_type(std::string);
        void set_path_output_spectra(std::string);
        void set_generate_report(int);
        void set_path_report(std::string);
        void set_generate_figure(int);
        void set_path_figure(std::string);
        void set_path_output_trend(std::string);
        void set_derivatives(int);
        void set_path_measurements(std::string);
        void set_path_input_spectrum(std::string);
        void set_path_energy_bins(std::string);
        void set_path_system_response(std::string);
        void set_path_icrp_factors(std::string);
        void set_path_ref_spectrum(std::string);
};


// class MeasurementManager {
//     public:
//         int num_meas_per_sample;

//         std::vector<double> measurements;
//         std::vector<double> measurements_scaled;
//         std::vector<double> std_devs;
// }

//--------------------------------------------------------------------------------------------------
// This class is designed to handle a single upper or lower uncertainty of a neutron fluence
// spectrum and the corresponding ambient dose equivalent. Implemented in a class instead of 
// directly in the unfolding method to allow access to variables later, for example when generating
// the report.
// This specific uncertainty manager class (J) is designed to calculate uncertainties using
// upper & lower boundaries on J, as established via a user-defined sigma value on J.
//--------------------------------------------------------------------------------------------------
class UncertaintyManagerJ {
    public:
        double j_threshold;
        double j_factor;
        int num_iterations;
        double dose_uncertainty;
        std::vector<double> bound_spectrum;
        std::vector<double> spectrum_uncertainty;

        UncertaintyManagerJ();
        UncertaintyManagerJ(double original_j_threshold, double sigma_j);

        void determineSpectrumUncertainty(std::vector<double> &mlemstop_spectrum, 
            int cutoff, int num_measurements, int num_bins, std::vector<double> &measurements, 
            std::vector<std::vector<double>> &nns_response, std::vector<double> &normalized_response,
            std::vector<double> &initial_spectrum);

        void determineDoseUncertainty(double dose, std::vector<double> &mlemstop_spectrum, int num_bins, 
            std::vector<double> &icrp_factors);
};


class UnfoldingReport {
    public:
        const std::string HEADER_DIVIDE = 
            "************************************************************************************************************************\n";
        const std::string SECTION_DIVIDE = 
            "\n========================================================================================================================\n\n";
        const std::string COLSTRING = 
            "--------------------";

        const int sw = 30; // settings column width
        const int cw = 20; // data column width
        const int rw = 9; // NNS response column width

        std::string path;
        std::string irradiation_conditions;
        std::vector<std::string> input_files;
        std::vector<std::string> input_file_flags;

        int cutoff;
        double error;
        double norm;
        double f_factor;
        int num_measurements;
        int num_bins;
        std::string uncertainty_type;
        int num_uncertainty_samples;
        std::string git_commit;

        std::vector<double> measurements; // measurements_report
        std::vector<double> measurements_nc;
        double dose_mu;
        double doserate_mu;
        int duration;
        std::string meas_units;

        std::vector<double> initial_spectrum;
        std::vector<double> energy_bins;
        std::vector<std::vector<double>> nns_response;
        std::vector<double> icrp_factors;

        std::vector<double> spectrum;
        std::vector<double> spectrum_uncertainty_upper;
        std::vector<double> spectrum_uncertainty_lower;
        int num_iterations;
        std::vector<double> mlem_ratio;
        double dose;
        double dose_uncertainty_upper;
        double dose_uncertainty_lower;
        double total_flux;
        double total_flux_uncertainty_upper;
        double total_flux_uncertainty_lower;
        double avg_energy;
        double avg_energy_uncertainty_upper;
        double avg_energy_uncertainty_lower;

        //MAP
        double beta;
        std::string algorithm;

        // MLEM-STOP
        int cps_crossover;
        double j_threshold;
        double j_final;
        UncertaintyManagerJ j_manager_low;
        UncertaintyManagerJ j_manager_high;
        int num_toss;

        UnfoldingReport(); 

        void prepare_report();
        void report_header(std::ofstream&);
        void report_settings(std::ofstream&);
        void report_measurement_info(std::ofstream&);
        void report_inputs(std::ofstream&);
        void report_mlem_info(std::ofstream&);
        void report_results(std::ofstream&);

        void set_path(std::string);
        void set_irradiation_conditions(std::string);
        void set_input_files(std::vector<std::string>&);
        void set_input_file_flags(std::vector<std::string>&);

        void set_cutoff(int);
        void set_error(double);
        void set_norm(double);
        void set_f_factor(double);
        void set_num_measurements(int);
        void set_num_bins(int);
        void set_uncertainty_type(std::string);
        void set_num_uncertainty_samples(int);
        void set_git_commit(std::string);

        void set_measurements(std::vector<double>&);
        void set_measurements_nc(std::vector<double>&);
        void set_dose_mu(double);
        void set_doserate_mu(double);
        void set_duration(int);
        void set_meas_units(std::string);

        void set_initial_spectrum(std::vector<double>&);
        void set_energy_bins(std::vector<double>&);
        void set_nns_response(std::vector<std::vector<double>>&);
        void set_icrp_factors(std::vector<double>&);

        void set_spectrum(std::vector<double>&);
        void set_spectrum_uncertainty_upper(std::vector<double>&);
        void set_spectrum_uncertainty_lower(std::vector<double>&);
        void set_num_iterations(int);
        void set_mlem_ratio(std::vector<double>&);
        void set_dose(double);
        void set_dose_uncertainty_upper(double);
        void set_dose_uncertainty_lower(double);
        void set_total_flux(double);
        void set_total_flux_uncertainty_upper(double);
        void set_total_flux_uncertainty_lower(double);
        void set_avg_energy(double);
        void set_avg_energy_uncertainty_upper(double);
        void set_avg_energy_uncertainty_lower(double);

        void set_algorithm(std::string);

        void set_cps_crossover(int);
        void set_j_threshold(double);
        void set_j_final(double);
        void set_j_manager_low(UncertaintyManagerJ);
        void set_j_manager_high(UncertaintyManagerJ);
        void set_num_toss(int);
};

class SpectraSettings{
    public:
        std::string path_input_data;
        std::string path_output_figure;
        std::string title;
        std::string x_label;
        std::string y_label;
        double x_min;
        double x_max;
        double y_min;
        double y_max;
        int x_res;
        int y_res;
        int y_num_divs;
        int y_digits_max;
        std::vector<std::string> legend_entries;
        std::vector<std::string> color_series;
        std::vector<std::string> color_error;
        int grayscale;
        std::vector<int> show_error;
        std::string error_style;
        int error_fill_style;
        int rows_per_spectrum;
        std::vector<int> line_style;
        std::vector<int> line_width;
        int border_width;
        int legend;
        std::vector<float> legend_coords;
        int textbox;
        std::vector<float> textbox_coords;
        std::vector<std::string> textbox_text;
        int plot_per_mu;
        std::vector<int> number_mu;
        std::vector<int> duration;
        int normalize;
        double margin_left;
        double margin_right;
        double margin_top;
        double margin_bottom;
        double font_size;
        double font_size_legend;
        double font_size_axis_labels;
        double font_size_axis_tick_labels;
        double font_size_title;
        double font_size_textbox;
        double x_label_offset;
        double y_label_offset;

        SpectraSettings();

        void set_setting(std::string,std::string);

        void set_path_input_data(std::string);
        void set_path_output_figure(std::string);
        void set_title(std::string);
        void set_x_label(std::string);
        void set_y_label(std::string);
        void set_x_min(std::string);
        void set_x_max(std::string);
        void set_y_min(std::string);
        void set_y_max(std::string);
        void set_x_res(std::string);
        void set_y_res(std::string);
        void set_y_num_divs(std::string);
        void set_y_digits_max(std::string);
        void set_legend_entries(std::string);
        void set_color_series(std::string);
        void set_color_error(std::string);
        void set_grayscale(std::string);
        void set_show_error(std::string);
        void set_error_style(std::string);
        void set_error_fill_style(std::string);
        void set_rows_per_spectrum(std::string);
        void set_line_style(std::string);
        void set_line_width(std::string);
        void set_border_width(std::string);
        void set_legend(std::string);
        void set_legend_coords(std::string);
        void set_textbox(std::string);
        void set_textbox_coords(std::string);
        void set_textbox_text(std::string);
        void set_plot_per_mu(std::string);
        void set_number_mu(std::string);
        void set_duration(std::string);
        void set_normalize(std::string);
        void set_margin_left(std::string);
        void set_margin_right(std::string);
        void set_margin_top(std::string);
        void set_margin_bottom(std::string);
        void set_font_size(std::string);
        void set_font_size_legend(std::string);
        void set_font_size_axis_labels(std::string);
        void set_font_size_axis_tick_labels(std::string);
        void set_font_size_title(std::string);
        void set_font_size_textbox(std::string);
        void set_x_label_offset(std::string);
        void set_y_label_offset(std::string);
};

class PlotSettings{
    public:
        std::string path_input_data;
        std::string path_output_figure;
        std::string data_format;
        std::string title;
        std::string x_label;
        std::string y_label;
        double y_min;
        double y_max;
        double x_min;
        double x_max;
        int x_res;
        int y_res;
        int x_num_divs;
        int y_num_divs;
        int x_log;
        int y_log;
        int x_grid;
        int y_grid;
        std::vector<std::string> legend_entries;
        std::vector<std::string> color_series;
        int grayscale;
        // std::vector<std::string> color_error;
        // std::vector<int> show_error;
        std::vector<int> line_style;
        std::vector<int> line_width;
        int border_width;
        int legend;
        std::vector<float> legend_coords;
        int textbox;
        std::vector<float> textbox_coords;
        std::vector<std::string> textbox_text;

        std::vector<std::string> plot_type;
        int legend_border_size;
        double legend_text_size;
        std::vector<int> marker_style;
        std::vector<int> marker_size;
        double margin_left;
        double margin_right;
        double margin_top;
        double margin_bottom;
        double x_label_offset;
        double y_label_offset;

        PlotSettings(); 

        void set_setting(std::string,std::string);

        void set_path_input_data(std::string);
        void set_path_output_figure(std::string);
        void set_data_format(std::string);
        void set_title(std::string);
        void set_x_label(std::string);
        void set_y_label(std::string);
        void set_y_min(std::string);
        void set_y_max(std::string);
        void set_x_min(std::string);
        void set_x_max(std::string);
        void set_x_res(std::string);
        void set_y_res(std::string);
        void set_x_log(std::string);
        void set_y_log(std::string);
        void set_x_grid(std::string);
        void set_y_grid(std::string);
        void set_x_num_divs(std::string);
        void set_y_num_divs(std::string);
        void set_legend_entries(std::string);
        void set_color_series(std::string);
        void set_grayscale(std::string);
        // void set_color_error(std::string);
        // void set_show_error(std::string);
        void set_line_style(std::string);
        void set_line_width(std::string);
        void set_border_width(std::string);
        void set_legend(std::string);
        void set_legend_coords(std::string);
        void set_textbox(std::string);
        void set_textbox_coords(std::string);
        void set_textbox_text(std::string);

        void set_plot_type(std::string);
        void set_legend_border_size(std::string);
        void set_legend_text_size(std::string);
        void set_marker_style(std::string);
        void set_marker_size(std::string);
        void set_margin_left(std::string);
        void set_margin_right(std::string);
        void set_margin_top(std::string);
        void set_margin_bottom(std::string);
        void set_x_label_offset(std::string);
        void set_y_label_offset(std::string);
};

class SurfaceSettings{
    public:
        std::string path_input_data;
        std::string path_output_figure;
        std::string title;
        std::string x_label;
        std::string y_label;
        std::string z_label;
        double y_min;
        double y_max;
        int x_min;
        int x_max;
        double z_min;
        double z_max;
        int x_res;
        int y_res;
        int x_num_divs;
        int z_num_divs;
        int color_palette;
        int num_color_bins;
        int border_width;

        SurfaceSettings(); 

        void set_setting(std::string,std::string);

        void set_path_input_data(std::string);
        void set_path_output_figure(std::string);
        void set_title(std::string);
        void set_x_label(std::string);
        void set_y_label(std::string);
        void set_z_label(std::string);
        void set_y_min(std::string);
        void set_y_max(std::string);
        void set_x_min(std::string);
        void set_x_max(std::string);
        void set_z_min(std::string);
        void set_z_max(std::string);
        void set_x_res(std::string);
        void set_y_res(std::string);
        void set_x_num_divs(std::string);
        void set_z_num_divs(std::string);
        void set_color_palette(std::string);
        void set_num_color_bins(std::string);
        void set_border_width(std::string);
};

#endif