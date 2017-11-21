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

#include "fileio.h"
#include "root_helpers.h"

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
#include "TVector.h"
#include "TLatex.h"

// Prototypes of helper functions
void stringToSVector(std::string test_string, std::vector<std::string>& result_vector);
void stringToIVector(std::string test_string, std::vector<int>& result_vector);

// This map variable contains all configurable settings for generating the output plot
// Default values are provided where approriate, but each of these settings may be configured
// by the user via the input/plot.cfg file
std::map<std::string,std::string> settings = {
    {"input_filename", "output_spectra.csv"},
    {"input_dir", "output/"},
    {"output_filename", "neutron_spectra.png"},
    {"output_dir", "output/"},
    {"title", "Neutron Flux Spectra"},
    {"x_label", "Energy / MeV"},
    {"y_label", "Fluence Rate / n#upoint cm^{-2}s^{-1}"},
    {"y_min", ""},
    {"y_max", ""},
    // {"x_tick_maj", ""},
    // {"x_tick_min", ""},
    // {"y_tick_maj", ""},
    // {"x_tick_min", ""},
    {"x_res", "800"},
    {"y_res", "600"},
    {"legend_entries", ""},
    {"color_spectra", "1,860,632,416"},
    {"color_error", "920,856,623,409"},
    {"show_error", "1,1,1,1"}, // Display uncertainty or not; 1=yes, 0=no
};

int main(int argc, char* argv[])
{
    // Set Settings
    std::string settings_file = "input/plot.cfg";
    setPlotSettings(settings_file, settings);
    std::string input_file = settings["input_dir"]+settings["input_filename"];
    std::string output_file = settings["output_dir"]+settings["output_filename"];

    // Read in data
    std::vector<std::string> headers;
    std::vector<double> energy_bins;
    std::vector<std::vector<double>> spectra_array;
    std::vector<std::vector<double>> error_array;

    readSpectra(input_file, headers, energy_bins, spectra_array, error_array);

    int num_spectra = spectra_array.size();
    int num_bins = energy_bins.size() - 1;

    // Convert from vector to array (array is necessary in TH1F constructors)
    double bins_arr[num_bins];
    for (int i_bin = 0; i_bin <= num_bins; i_bin++)
    {
        bins_arr[i_bin] = energy_bins[i_bin];
    }

    // Generate the plot
    TCanvas *c1 = new TCanvas("c1","c1",atoi(settings["x_res"].c_str()),atoi(settings["y_res"].c_str())); // Resulution of the graph (px) specified in parameters
    TLegend* leg = new TLegend(0.15, 0.7, 0.48, 0.85); // startx, starty, endx, endy

    // Convert string settings (comma-delimited) to vectors
    std::vector<std::string> legend_entries;
    std::vector<int> color_spectra;
    std::vector<int> color_error;
    std::vector<int> show_error;
    if (!settings["legend_entries"].empty())
        stringToSVector(settings["legend_entries"], legend_entries);
    if (!settings["color_spectra"].empty())
        stringToIVector(settings["color_spectra"], color_spectra);
    if (!settings["color_error"].empty())
        stringToIVector(settings["color_error"], color_error);
    if (!settings["show_error"].empty())
        stringToIVector(settings["show_error"], show_error);

    // Create histograms
    std::vector<TH1F*> histograms;
    for (int i_spec = 0; i_spec < num_spectra; i_spec++) {
        histograms.push_back(new TH1F("","",num_bins,bins_arr));
    }
    for (int i_spec = 0; i_spec < num_spectra; i_spec++) {
        // Fill the histograms with spectral data
        for (int i_bin = 0; i_bin < num_bins; i_bin++) {
            histograms[i_spec]->Fill(energy_bins[i_bin], spectra_array[i_spec][i_bin]);
        }

        // Apply settings
        histograms[i_spec]->SetStats(0);   // Do not show the stats (mean and standard deviation);
        histograms[i_spec]->SetLineWidth(1);
        histograms[i_spec]->SetLineColor(color_spectra[i_spec%color_spectra.size()]);

        // Add current spectrum to legend
        if (legend_entries.empty())
            leg->AddEntry(histograms[i_spec], headers[i_spec].c_str(), "l");
        else
            leg->AddEntry(histograms[i_spec], legend_entries[i_spec%legend_entries.size()].c_str(), "l");

        // Title & axes manipulations for first spectrum
        if (i_spec == 0) {
            histograms[i_spec]->SetTitle(settings["title"].c_str());
            histograms[i_spec]->GetXaxis()->SetTitleOffset(1.4);
            histograms[i_spec]->GetXaxis()->CenterTitle();
            histograms[i_spec]->SetXTitle(settings["x_label"].c_str());
            histograms[i_spec]->GetYaxis()->SetTitleOffset(1.4);
            histograms[i_spec]->GetYaxis()->CenterTitle();
            histograms[i_spec]->SetYTitle(settings["y_label"].c_str());

            histograms[i_spec]->SetTickLength(0.015,"xy"); // Length of tick marks (default = 0.02)

            if (!settings["y_min"].empty())
                histograms[i_spec]->SetMinimum(atoi(settings["y_min"].c_str()));
            if (!settings["y_max"].empty())
                histograms[i_spec]->SetMaximum(atoi(settings["y_max"].c_str()));

            histograms[i_spec]->Draw("HIST");  // Draw the histogram without the error bars;
        }
        // Plot subsequent spectra on the same axes
        else {
            histograms[i_spec]->Draw("HIST SAME");
        }
    }

    leg->Draw(); // Add the legend to the canvas

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

//==================================================================================================
// Convert a comma-delimited string into a vector of strings.
//
// Args:
//  - test_string: the comma-delimited string to be processed
//  - result_vector: the vector that will be assigned string values
//==================================================================================================
void stringToSVector(std::string test_string, std::vector<std::string>& result_vector) {
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
    std::istringstream line_stream(test_string);
    std::string stoken; // store individual values between delimiters on a line

    // Loop through each line, delimiting at commas
    while (getline(line_stream, stoken, ',')) {
        result_vector.push_back(atoi(stoken.c_str()));
    }
}