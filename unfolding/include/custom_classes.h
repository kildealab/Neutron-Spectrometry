#ifndef CUSTOM_CLASSES_H
#define CUSTOM_CLASSES_H

#include <stdlib.h>
#include <string>
#include <vector>

class UnfoldingSettings {
    public:
        double norm; // vendor specfied normalization factor for the NNS used
        double error; // The target error on ratio between the experimental data points and the estimated data points from unfolding (e.g. 0.1 means the values must be within 10% of each other before the algorithm will terminate)
        double f_factor; // factor that converts measured charge (in fC) to counts per second [fA/cps]
        int cutoff; // maximum # of iterations in the unfolding algorithm
        int num_poisson_samples; // # of samples from Poisson distribution for uncertainty estimation
        // MAP specific
        double beta; // factor that affects damping of high noise spectral component in MAP method
        std::string prior;
        // Optimize specific
        int min_num_iterations;
        int max_num_iterations;
        int iteration_increment;
        double min_beta;
        double max_beta;
        std::string parameter_of_interest;
        std::string algorithm;
        std::string trend_type;
        std::string auto_output_path;

        UnfoldingSettings(); 

        void set_setting(std::string,std::string);

        void set_norm(double);
        void set_error(double);
        void set_f_factor(double);
        void set_cutoff(int);
        void set_num_poisson_samples(int);
        void set_beta(double);
        void set_prior(std::string);
        void set_min_num_iterations(int);
        void set_max_num_iterations(int);
        void set_iteration_increment(int);
        void set_min_beta(double);
        void set_max_beta(double);
        void set_parameter_of_interest(std::string);
        void set_algorithm(std::string);
        void set_trend_type(std::string);
        void set_auto_output_path(std::string);
};

class SpectraSettings {
    public:
        std::string input_filename;
        std::string input_dir;
        std::string output_filename;
        std::string output_dir;
        std::string title;
        std::string x_label;
        std::string y_label;
        double y_min;
        double y_max;
        int x_res;
        int y_res;
        int y_num_divs;
        std::vector<std::string> legend_entries;
        std::vector<std::string> color_series;
        std::vector<std::string> color_error;
        std::vector<int> show_error;
        std::vector<int> line_style;
        std::vector<int> line_width;
        std::vector<float> legend_coords;
        int textbox;
        std::vector<float> textbox_coords;
        std::vector<std::string> textbox_text;
        int plot_per_mu;
        std::vector<int> number_mu;
        std::vector<int> duration;

        SpectraSettings();

        void set_setting(std::string,std::string);

        void set_input_filename(std::string);
        void set_input_dir(std::string);
        void set_output_filename(std::string);
        void set_output_dir(std::string);
        void set_title(std::string);
        void set_x_label(std::string);
        void set_y_label(std::string);
        void set_y_min(std::string);
        void set_y_max(std::string);
        void set_x_res(std::string);
        void set_y_res(std::string);
        void set_y_num_divs(std::string);
        void set_legend_entries(std::string);
        void set_color_series(std::string);
        void set_color_error(std::string);
        void set_show_error(std::string);
        void set_line_style(std::string);
        void set_line_width(std::string);
        void set_legend_coords(std::string);
        void set_textbox(std::string);
        void set_textbox_coords(std::string);
        void set_textbox_text(std::string);
        void set_plot_per_mu(std::string);
        void set_number_mu(std::string);
        void set_duration(std::string);
};

class PlotSettings {
    public:
        std::string input_filename;
        std::string input_dir;
        std::string output_filename;
        std::string output_dir;
        std::string title;
        std::string x_label;
        std::string y_label;
        double y_min;
        double y_max;
        int x_min;
        int x_max;
        int x_res;
        int y_res;
        int x_num_divs;
        int y_num_divs;
        std::vector<std::string> legend_entries;
        std::vector<std::string> color_series;
        // std::vector<std::string> color_error;
        // std::vector<int> show_error;
        std::vector<int> line_style;
        std::vector<int> line_width;
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

        void set_input_filename(std::string);
        void set_input_dir(std::string);
        void set_output_filename(std::string);
        void set_output_dir(std::string);
        void set_title(std::string);
        void set_x_label(std::string);
        void set_y_label(std::string);
        void set_y_min(std::string);
        void set_y_max(std::string);
        void set_x_min(std::string);
        void set_x_max(std::string);
        void set_x_res(std::string);
        void set_y_res(std::string);
        void set_x_num_divs(std::string);
        void set_y_num_divs(std::string);
        void set_legend_entries(std::string);
        void set_color_series(std::string);
        // void set_color_error(std::string);
        // void set_show_error(std::string);
        void set_line_style(std::string);
        void set_line_width(std::string);
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

class SurfaceSettings {
    public:
        std::string input_filename;
        std::string input_dir;
        std::string output_filename;
        std::string output_dir;
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

        SurfaceSettings(); 

        void set_setting(std::string,std::string);

        void set_input_filename(std::string);
        void set_input_dir(std::string);
        void set_output_filename(std::string);
        void set_output_dir(std::string);
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
};

#endif