//**************************************************************************************************
// This program reads in an arbitrary number of spectra and their uncertainties from a text file and
// plots them on a single set of axes.
//**************************************************************************************************

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <map>

// #include "custom_classes.h"
#include "fileio.h"
#include "root_helpers.h"

// Root
#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TFrame.h"
#include "TGaxis.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TText.h"
#include "TVector.h"
#include "TVirtualPad.h"

// Prototypes of helper functions
// void stringToSVector(std::string test_string, std::vector<std::string>& result_vector);
// void stringToIVector(std::string test_string, std::vector<int>& result_vector);
// void stringToFVector(std::string test_string, std::vector<float>& result_vector);

// This map variable contains all configurable settings for generating the output plot
// Default values are provided where approriate, but each of these settings may be configured
// by the user via the input/plot.cfg file
std::map<std::string,std::string> settings = {
    {"input_filename", "output_spectra.csv"},
    {"input_dir", "output/"},
    {"output_filename", "neutron_spectra.png"},
    {"output_dir", "output/"},
    {"title", "Neutron Fluence Spectra"},
    {"x_label", "Energy (MeV)"},
    {"y_label", "Fluence (n #upoint cm^{-2} s^{-1})"},
    // {"y_label", "Fluence / n#upoint cm^{-2} MU^{-1}"},
    {"y_min", ""},
    {"y_max", ""},
    // {"x_tick_maj", ""},
    // {"x_tick_min", ""},
    // {"y_tick_maj", ""},
    // {"x_tick_min", ""},
    {"x_res", "800"},
    {"y_res", "600"},
    // For instructions on number of divisions, refer to https://root.cern.ch/doc/master/classTAttAxis.html
    {"y_num_divs", ""}, 
    {"legend_entries", ""},
    {"color_spectra", "#000000,#C63822,#607FD5,#55A961"}, // black,red,blue,green
    {"color_error", "#333333,#E79A9F,#6B8EF0,#69CF77"}, // slightly different shades
    {"show_error", "1,1,1,1"}, // Display uncertainty or not; 1=yes, 0=no
    // For line styling, refer to https://root.cern.ch/doc/master/classTAttLine.html
    {"line_style", "1,1,1,1"},
    {"line_width", "1,1,1,1"},
    {"legend_coords", "0.15,0.75,0.4,0.85"}, // start_x, start_y, end_x, end_y
    {"textbox", "0"}, // 0 = false, 1 = true
    {"textbox_coords", "0,0,0,0"}, // start_x, start_y, end_x, end_y
    {"textbox_text", ""}, // Each line of text separated by commas
    {"plot_per_mu", "0"}, // 0 = false, 1 = true
    {"number_mu", ""},
    {"duration", ""},
};

int main(int argc, char* argv[])
{
    // Set Settings
    std::string settings_file = "input/plot.cfg";
    setPlotSettingsOld(settings_file, settings); // Fill settings map with any user provided settings
    std::string input_file = settings["input_dir"]+settings["input_filename"];
    std::string output_file = settings["output_dir"]+settings["output_filename"];

    // Convert string settings (comma-delimited) to vectors of appropriate type
    std::vector<std::string> legend_entries;
    std::vector<std::string> color_spectra_temp; // Temporarily hold string that represents Hex color
    std::vector<std::string> color_error_temp; // Temporarily hold string that represents Hex color
    std::vector<int> color_spectra;
    std::vector<int> color_error;
    std::vector<int> show_error;
    std::vector<int> line_style;
    std::vector<int> line_width;
    std::vector<float> legend_coords;
    bool textbox = atoi(settings["textbox"].c_str());
    std::vector<float> textbox_coords;
    std::vector<std::string> textbox_text;
    bool plot_per_mu = atoi(settings["plot_per_mu"].c_str());
    std::vector<int> number_mu;
    std::vector<int> duration;

    if (!settings["legend_entries"].empty())
        stringToSVector(settings["legend_entries"], legend_entries);
    if (!settings["show_error"].empty())
        stringToIVector(settings["show_error"], show_error);
    if (!settings["line_style"].empty())
        stringToIVector(settings["line_style"], line_style);
    if (!settings["line_width"].empty())
        stringToIVector(settings["line_width"], line_width);
    if (!settings["legend_coords"].empty())
        stringToDVector(settings["legend_coords"], legend_coords);
    if (!settings["textbox_coords"].empty())
        stringToDVector(settings["textbox_coords"], textbox_coords);
    if (!settings["textbox_text"].empty())
        stringToSVector(settings["textbox_text"], textbox_text);
    if (!settings["number_mu"].empty())
        stringToIVector(settings["number_mu"], number_mu);
    if (!settings["duration"].empty())
        stringToIVector(settings["duration"], duration);

    // The color settings must be converted from strings to integers corresponding to the provided Hex color
    // The TColor:GetColor method converts from hex to TColor (which is an int)
    if (!settings["color_spectra"].empty()) {
        stringToSVector(settings["color_spectra"], color_spectra_temp); // Break apart string into vector of strings
        for (int i=0; i < color_spectra_temp.size(); i++){
            color_spectra.push_back(TColor::GetColor(color_spectra_temp[i].c_str())); // Convert hex color (stored as string) to TColor color (stored as int)
        }
    }
    if (!settings["color_error"].empty()) {
        stringToSVector(settings["color_error"], color_error_temp); // Break apart string into vector of strings
        for (int i=0; i < color_error_temp.size(); i++){
            color_error.push_back(TColor::GetColor(color_error_temp[i].c_str())); // Convert hex color (stored as string) to TColor color (stored as int)
        }
    }

    // Adjust y-axis label if plotting per MU
    if (plot_per_mu){
        settings["y_label"] = "Fluence (n #upoint cm^{-2} MU^{-1})";
    }

    // Read in data
    std::vector<std::string> headers;
    std::vector<double> energy_bins;
    std::vector<std::vector<double>> spectra_array;
    std::vector<std::vector<double>> error_array;
    readSpectra(input_file, headers, energy_bins, spectra_array, error_array, plot_per_mu, number_mu, duration);

    int num_spectra = spectra_array.size();
    int num_bins = energy_bins.size() - 1;


    // Generate the plot area
    TCanvas *c1 = new TCanvas("c1","c1",atoi(settings["x_res"].c_str()),atoi(settings["y_res"].c_str())); // Resulution of the graph (px) specified in parameters

    // Generate the legend
    TLegend* leg = new TLegend(legend_coords[0], legend_coords[1], legend_coords[2], legend_coords[3]); // with a text box
    leg->SetBorderSize(0);
    leg->SetTextSize(0.035);

    // Convert spectral data from vector to array (array is necessary in TH1F constructors)
    double bins_arr[num_bins];
    for (int i_bin = 0; i_bin <= num_bins; i_bin++)
    {
        bins_arr[i_bin] = energy_bins[i_bin];
    }

    // Create histograms
    std::vector<TH1F*> histograms;
    for (int i_spec = 0; i_spec < num_spectra; i_spec++) {
        histograms.push_back(new TH1F("","",num_bins,bins_arr));
    }

    // Apply settings to the newly-created histograms
    for (int i_spec = 0; i_spec < num_spectra; i_spec++) {
        // Fill the histograms with spectral data
        for (int i_bin = 0; i_bin < num_bins; i_bin++) {
            histograms[i_spec]->Fill(energy_bins[i_bin], spectra_array[i_spec][i_bin]);
        }

        // Apply settings
        histograms[i_spec]->SetStats(0);   // Do not show the stats (mean and standard deviation);
        histograms[i_spec]->SetLineColor(color_spectra[i_spec%color_spectra.size()]);
        histograms[i_spec]->SetLineStyle(line_style[i_spec%line_style.size()]);
        histograms[i_spec]->SetLineWidth(line_width[i_spec%line_width.size()]);

        // Add current spectrum to legend
        // If user did not provide legend titles, use the "header" names provided in the input spectrum file
        if (legend_entries.empty())
            leg->AddEntry(histograms[i_spec], headers[i_spec].c_str(), "l");
        else
            leg->AddEntry(histograms[i_spec], legend_entries[i_spec%legend_entries.size()].c_str(), "l");
        // if (i_spec == 1) {
        //     leg->AddEntry((TObject*)0, "", "");
        // }

        // Title & axes manipulations. Only needs to be done for first spectrum
        if (i_spec == 0) {
            histograms[i_spec]->SetTitle(settings["title"].c_str());

            histograms[i_spec]->GetXaxis()->SetTitleOffset(1.4);
            histograms[i_spec]->GetXaxis()->CenterTitle();
            histograms[i_spec]->SetXTitle(settings["x_label"].c_str());

            histograms[i_spec]->GetYaxis()->SetTitleOffset(1.4);
            histograms[i_spec]->GetYaxis()->CenterTitle();
            histograms[i_spec]->SetYTitle(settings["y_label"].c_str());

            histograms[i_spec]->SetTickLength(0.015,"xy"); // Length of tick marks (default = 0.02)

            // Set x-axis range
            histograms[i_spec]->GetXaxis()->SetRange(0,52); // Range according to bin number
            // histograms[i_spec]->GetXaxis()->SetRangeUser(1e-9,10); // Range according to value (must be in range spanned by bins)
            
            // Set axis tick label size
            histograms[i_spec]->GetYaxis()->SetLabelSize(0.03);
            histograms[i_spec]->GetXaxis()->SetLabelSize(0.03);

            if (!settings["y_min"].empty())
                histograms[i_spec]->SetMinimum(atof(settings["y_min"].c_str()));
            if (!settings["y_max"].empty())
                histograms[i_spec]->SetMaximum(atof(settings["y_max"].c_str()));
            // Set # of divisions on the y axis
            // refer to https://root.cern.ch/doc/master/classTAttAxis.html#ae3067b6d4218970d09418291cbd84084
            if (!settings["y_num_divs"].empty())
                histograms[i_spec]->GetYaxis()->SetNdivisions(atoi(settings["y_num_divs"].c_str()),kFALSE); // kFALSE forces ROOT to use the number of divisions specified

            histograms[i_spec]->Draw("HIST");  // Draw the histogram without the error bars;
        }
        // Plot subsequent spectra on the same axes
        else {
            histograms[i_spec]->Draw("HIST SAME");
        }
    }

    // Add the legend to the canvas
    leg->Draw();

    // Optionally add a text box. Must draw the text box later in script as well, if necessary
    if(textbox){
        // TPaveText* pt = new TPaveText(0.15, 0.75, 0.4, 0.85, "nbNDC"); // nb specifies no border, NDC specifies method of defining coordinates
        TPaveText* pt = new TPaveText(textbox_coords[0], textbox_coords[1], textbox_coords[2], textbox_coords[3], "nbNDC"); // nb specifies no border, NDC specifies method of defining coordinates
        pt->SetFillColorAlpha(kWhite,1);
        pt->SetTextAlign(12);
        pt->SetTextSize(0.035);
        pt->SetTextFont(42);

        for (int i=0; i<textbox_text.size(); i++){
            pt->AddText(textbox_text[i].c_str());
        }
        
        pt->Draw(); // Add the text box
    }

    // Prepare the Uncertainties
    // Need middle energy bin values at which the uncertainties are plotted 
    std::vector<double> bins_avg;
    for (int i_bin = 0; i_bin < num_bins; i_bin++)
    {
        bins_avg.push_back((energy_bins[i_bin] + energy_bins[i_bin+1])/2);
    }

    // Plot the uncertainties
    for (int i_spec = 0; i_spec < num_spectra; i_spec++) {
        TGraphErrors *ge = new TGraphErrors(num_bins, &(bins_avg[0]), &(spectra_array[i_spec][0]), 0, &(error_array[i_spec][0]));
        ge->SetFillColor(color_error[i_spec%color_error.size()]);
        ge->SetFillStyle(3002); // code for fill pattern

        // Add the uncertainty to the canvas if set as such
        if (show_error[i_spec%show_error.size()])
            ge->Draw("P3");
    }

    // Prepare the canvas for output
    TGaxis::SetMaxDigits(3); // Set maximum number of digits (i.e. for scientific notation) on yaxis
    c1->SetLogx();
    c1->Update();
    c1->Modified();
    c1->SetTickx(); // No parameter means show tick marks on both sides, labels on one
    c1->SetTicky(); // No parameter means show tick marks on both sides, labels on one

    // Output the plot to file
    const char *cstr_figure_file = output_file.c_str();
    c1->Print(cstr_figure_file);

    return 0;
}