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

#include "custom_classes.h"
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

int main(int argc, char* argv[])
{
    // Set Settings
    std::string settings_file = "input/plot_spectra.cfg";
    SpectraSettings settings;
    setSpectraSettings(settings_file, settings);

    std::string input_file = settings.input_dir + settings.input_filename;
    std::string output_file = settings.output_dir + settings.output_filename;

    // Read in data
    std::vector<std::string> headers;
    std::vector<double> energy_bins;
    std::vector<std::vector<double>> spectra_array;
    std::vector<std::vector<double>> error_array;
    readSpectra(input_file, headers, energy_bins, spectra_array, error_array, settings.plot_per_mu, settings.number_mu, settings.duration);

    int num_spectra = spectra_array.size();
    int num_bins = energy_bins.size() - 1;

    // If spectra should be normalized
    if (settings.normalize) {
        for (int i_spec=0; i_spec < num_spectra; i_spec++) {
            // Get max value
            int max_value = 0;
            for (int i_bin=0; i_bin < num_bins; i_bin++) {
                if (spectra_array[i_spec][i_bin] > max_value) {
                    max_value = spectra_array[i_spec][i_bin];
                }
            }
            // Normalize
            for (int i_bin=0; i_bin < num_bins; i_bin++) {
                spectra_array[i_spec][i_bin] = spectra_array[i_spec][i_bin] / max_value;
                error_array[i_spec][i_bin] = error_array[i_spec][i_bin] / max_value;
            }
        }
    }

    // Generate the plot area
    TCanvas *c1 = new TCanvas("c1","c1",settings.x_res,settings.y_res); // Resolution of the graph (px) specified in parameters

    // Generate the legend
    TLegend* leg = new TLegend(settings.legend_coords[0], settings.legend_coords[1], settings.legend_coords[2], settings.legend_coords[3]); // with a text box
    leg->SetBorderSize(0);
    leg->SetTextSize(0.035);

    // Convert spectral data from vector to array (array is necessary in TH1F constructors)
    double bins_arr[num_bins];
    for (int i_bin = 0; i_bin <= num_bins; i_bin++)
    {
        bins_arr[i_bin] = energy_bins[i_bin];
    }

    int setting_size; // normalize all settings by the number of entries in that setting (e.g. if have 10 colors set, the 11th plot series should reuse the 1st color)


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
        setting_size = settings.color_series.size();
        histograms[i_spec]->SetLineColor(TColor::GetColor(settings.color_series[i_spec%setting_size].c_str()));

        // width
        setting_size = settings.line_width.size();
        histograms[i_spec]->SetLineWidth(settings.line_width[i_spec%setting_size]);
        // style
        setting_size = settings.line_style.size();
        histograms[i_spec]->SetLineStyle(settings.line_style[i_spec%setting_size]);


        // Add current spectrum to legend
        // If user did not provide legend titles, use the "header" names provided in the input spectrum file
        if (settings.legend_entries.empty()) {
            leg->AddEntry(histograms[i_spec], headers[i_spec].c_str(), "l");
        }
        else {
            leg->AddEntry(histograms[i_spec], settings.legend_entries[i_spec].c_str(), "l");
        }
        // if (i_spec == 1) {
        //     leg->AddEntry((TObject*)0, "", "");
        // }

        // Title & axes manipulations. Only needs to be done for first spectrum
        if (i_spec == 0) {
            histograms[i_spec]->SetTitle(settings.title.c_str());

            histograms[i_spec]->GetXaxis()->SetTitleOffset(1.4);
            histograms[i_spec]->GetXaxis()->CenterTitle();
            histograms[i_spec]->SetXTitle(settings.x_label.c_str());

            histograms[i_spec]->GetYaxis()->SetTitleOffset(1.4);
            histograms[i_spec]->GetYaxis()->CenterTitle();
            histograms[i_spec]->SetYTitle(settings.y_label.c_str());

            histograms[i_spec]->SetTickLength(0.015,"xy"); // Length of tick marks (default = 0.02)

            // Set x-axis range
            histograms[i_spec]->GetXaxis()->SetRange(0,52); // Range according to bin number
            // histograms[i_spec]->GetXaxis()->SetRangeUser(1e-9,10); // Range according to value (must be in range spanned by bins)
            
            // Set axis tick label size
            histograms[i_spec]->GetYaxis()->SetLabelSize(0.03);
            histograms[i_spec]->GetXaxis()->SetLabelSize(0.03);

            if (settings.y_min != settings.y_max) {
                histograms[i_spec]->SetMinimum(settings.y_min);
                histograms[i_spec]->SetMaximum(settings.y_max);
            }

            // Set # of divisions on the y axis
            // refer to https://root.cern.ch/doc/master/classTAttAxis.html#ae3067b6d4218970d09418291cbd84084
            if (settings.y_num_divs != 0)
                histograms[i_spec]->GetYaxis()->SetNdivisions(settings.y_num_divs,kFALSE);


            histograms[i_spec]->Draw("HIST");  // Draw the histogram without the error bars;
        }
        // Plot subsequent spectra on the same axes
        else {
            histograms[i_spec]->Draw("HIST SAME");
        }
    }

    // Add the legend to the canvas
    if (settings.legend) {
        leg->Draw();
    }

    // Optionally add a text box. Must draw the text box later in script as well, if necessary
    if(settings.textbox){
        // TPaveText* pt = new TPaveText(0.15, 0.75, 0.4, 0.85, "nbNDC"); // nb specifies no border, NDC specifies method of defining coordinates
        TPaveText* pt = new TPaveText(settings.textbox_coords[0], settings.textbox_coords[1], settings.textbox_coords[2], settings.textbox_coords[3], "nbNDC"); // nb specifies no border, NDC specifies method of defining coordinates
        pt->SetFillColorAlpha(kWhite,1);
        pt->SetTextAlign(12);
        pt->SetTextSize(0.035);
        pt->SetTextFont(42);

        for (int i=0; i<settings.textbox_text.size(); i++){
            pt->AddText(settings.textbox_text[i].c_str());
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
        setting_size = settings.color_error.size();
        ge->SetFillColor(TColor::GetColor(settings.color_error[i_spec%setting_size].c_str()));

        ge->SetFillStyle(3002); // code for fill pattern

        // Add the uncertainty to the canvas if set as such
        setting_size = settings.show_error.size();
        if (settings.show_error[i_spec%setting_size])
            ge->Draw("P3");
    }

    // Prepare the canvas for output
    TGaxis::SetMaxDigits(3); // Set maximum number of digits (i.e. for scientific notation) on yaxis
    c1->SetLogx();
    c1->Update();
    c1->Modified();
    c1->SetTickx(); // No parameter means show tick marks on both sides, labels on one
    c1->SetTicky(); // No parameter means show tick marks on both sides, labels on one
    gStyle->SetLineWidth(settings.border_width); // Set width of axis/border around the plot

    // Output the plot to file
    const char *cstr_figure_file = output_file.c_str();
    c1->Print(cstr_figure_file);

    return 0;
}