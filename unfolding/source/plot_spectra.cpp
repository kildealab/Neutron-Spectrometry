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
#include "TGraphAsymmErrors.h"
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

    // Read in data
    std::vector<std::string> headers;
    std::vector<double> energy_bins;
    std::vector<std::vector<double>> spectra_array;
    std::vector<std::vector<double>> error_upper_array;
    std::vector<std::vector<double>> error_lower_array;
    readSpectra(settings.path_input_data, headers, energy_bins, spectra_array, error_lower_array, error_upper_array,
        settings.plot_per_mu, settings.number_mu, settings.duration, settings.rows_per_spectrum);

    int num_spectra = spectra_array.size();
    int num_bins = energy_bins.size() - 1;

    // If spectra should be normalized
    if (settings.normalize) {
        for (int i_spec=0; i_spec < num_spectra; i_spec++) {
            // Get max value
            double max_value = 0;
            for (int i_bin=0; i_bin < num_bins; i_bin++) {
                if (spectra_array[i_spec][i_bin] > max_value) {
                    max_value = spectra_array[i_spec][i_bin];
                }
            }
            // Normalize
            for (int i_bin=0; i_bin < num_bins; i_bin++) {
                spectra_array[i_spec][i_bin] = spectra_array[i_spec][i_bin] / max_value;
                error_lower_array[i_spec][i_bin] = error_lower_array[i_spec][i_bin] / max_value;
                error_upper_array[i_spec][i_bin] = error_upper_array[i_spec][i_bin] / max_value;
            }
        }
    }

    // Generate the plot area
    TCanvas *c1 = new TCanvas("c1","c1",settings.x_res,settings.y_res); // Resolution of the graph (px) specified in parameters
    if (settings.grayscale){
        c1->GetCanvas()->SetGrayscale();
    }

    // Create a ghost histogram that defines the x range. There are no other options in ROOT to
    // suitably set the x range of a TH1 object. See:
    // https://root-forum.cern.ch/t/th1-xaxis-setlimits/2880/3
    // https://root-forum.cern.ch/t/the-good-old-x-axis-range-problem/33840/6

    // x axis limits
    if (settings.x_min == settings.x_max) {
        settings.x_min = energy_bins[0];
        settings.x_max = energy_bins[num_bins];
    }

    TH1F *ghosthist = new TH1F ("","",1,settings.x_min,settings.x_max) ;
    ghosthist->SetBinContent(1,0); // ghosthist has one bin that extends from min to max, value 0
    ghosthist->SetLineWidth(0); // don't show the ghosthist bin

    c1->SetLeftMargin(settings.margin_left); 
    c1->SetRightMargin(settings.margin_right); 
    c1->SetTopMargin(settings.margin_top); 
    c1->SetBottomMargin(settings.margin_bottom); 

    // Title & axes manipulations:
    ghosthist->SetTitle(settings.title.c_str());

    ghosthist->GetXaxis()->SetTitleOffset(settings.x_label_offset);
    ghosthist->GetXaxis()->CenterTitle();
    ghosthist->SetXTitle(settings.x_label.c_str());
    ghosthist->GetXaxis()->SetTitleSize(0.035*settings.font_size*settings.font_size_axis_labels);

    ghosthist->GetYaxis()->SetTitleOffset(settings.y_label_offset);
    ghosthist->GetYaxis()->CenterTitle();
    ghosthist->SetYTitle(settings.y_label.c_str());
    ghosthist->GetYaxis()->SetTitleSize(0.035*settings.font_size*settings.font_size_axis_labels);

    ghosthist->SetTickLength(0.015,"xy"); // Length of tick marks (default = 0.02)
    
    // Set axis tick label size
    ghosthist->GetYaxis()->SetLabelSize(0.03*settings.font_size*settings.font_size_axis_tick_labels);
    ghosthist->GetXaxis()->SetLabelSize(0.03*settings.font_size*settings.font_size_axis_tick_labels);

    // y axis limits
    if (settings.y_min != settings.y_max) {
        ghosthist->SetMinimum(settings.y_min);
        ghosthist->SetMaximum(settings.y_max);
    }
    // Find max spectral value if no limits provided by user
    else {
        double global_max_value = 0;
        // Get max spectra value
        for (int i_spec=0; i_spec < num_spectra; i_spec++) {
            // Get max value
            double max_value = 0;
            for (int i_bin=0; i_bin < num_bins; i_bin++) {
                if (spectra_array[i_spec][i_bin] > max_value) {
                    max_value = spectra_array[i_spec][i_bin];
                }
            }
            if (max_value > global_max_value) {
                global_max_value = max_value;
            }
        }
        ghosthist->SetMaximum(global_max_value/0.8);
    }

    // Set # of divisions on the y axis
    // refer to https://root.cern.ch/doc/master/classTAttAxis.html#ae3067b6d4218970d09418291cbd84084
    if (settings.y_num_divs != 0)
        ghosthist->GetYaxis()->SetNdivisions(settings.y_num_divs,kFALSE);

    // Add ghost to the canvas
    ghosthist->Draw("HIST");
    // Adjust title size: for some reason it works here but not after the next line
    gStyle->SetTitleSize(0.05*settings.font_size*settings.font_size_title,"t");
    gStyle->SetOptStat(0); // Turns off the statistics box that is generated by default

    // Generate the legend
    TLegend* leg = new TLegend(
        settings.legend_coords[0], settings.legend_coords[1], settings.legend_coords[2], settings.legend_coords[3]
    ); // with a text box
    leg->SetBorderSize(0);
    leg->SetTextSize(0.035*settings.font_size*settings.font_size_legend);

    // Convert spectral data from vector to array (array is necessary in TH1F constructors)
    double bins_arr[num_bins];
    for (int i_bin = 0; i_bin <= num_bins; i_bin++)
    {
        bins_arr[i_bin] = energy_bins[i_bin];
    }

    int setting_size; // normalize all settings by the number of entries in that setting 
                      // (e.g. if have 10 colors set, the 11th plot series should reuse the 1st color)

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
        //     leg->AddEntry((TObject*)0, "", "");
        // }

    }

    // Add the legend to the canvas
    if (settings.legend) {
        leg->Draw();
    }

    // Optionally add a text box. Must draw the text box later in script as well, if necessary
    if(settings.textbox){
        // TPaveText* pt = new TPaveText(0.15, 0.75, 0.4, 0.85, "nbNDC"); 
        // nb specifies no border, NDC specifies method of defining coordinates
        TPaveText* pt = new TPaveText(settings.textbox_coords[0], settings.textbox_coords[1], 
            settings.textbox_coords[2], settings.textbox_coords[3], "nbNDC"
        );
        pt->SetFillColorAlpha(kWhite,0);
        pt->SetTextAlign(12);
        pt->SetTextSize(0.035*settings.font_size*settings.font_size_textbox);
        pt->SetTextFont(62);

        for (int i=0; i<settings.textbox_text.size(); i++){
            pt->AddText(settings.textbox_text[i].c_str());
        }
        
        pt->Draw(); // Add the text box
    }

    // Prepare the Uncertainties
    // Need middle energy bin values at which the uncertainties are plotted 
    std::vector<double> bins_avg;
    std::vector<double> bins_log_avg;
    std::vector<double> xerror_lower_array;
    std::vector<double> xerror_upper_array;
    for (int i_bin = 0; i_bin < num_bins; i_bin++)
    {
        bins_avg.push_back((energy_bins[i_bin] + energy_bins[i_bin+1])/2);
        double log_avg = sqrt(energy_bins[i_bin]*energy_bins[i_bin+1]);
        bins_log_avg.push_back(log_avg);
        xerror_lower_array.push_back(log_avg - energy_bins[i_bin]);
        xerror_upper_array.push_back(energy_bins[i_bin+1] - log_avg);
    }

    // Plot the uncertainties
    // Draw in inverse order so the first plot is "at the front" of the graph.
    for (int i_spec = num_spectra-1; i_spec >= 0; i_spec--) {
        // TGraphErrors *ge = new TGraphErrors(
        //     num_bins, &(bins_avg[0]), &(spectra_array[i_spec][0]), 0, &(error_array[i_spec][0])
        // );
        TGraphAsymmErrors *ge = new TGraphAsymmErrors(
            num_bins, &(bins_log_avg[0]), &(spectra_array[i_spec][0]), &(xerror_lower_array[0]), &(xerror_upper_array[0]), 
            &(error_lower_array[i_spec][0]), &(error_upper_array[i_spec][0])
        );

        setting_size = settings.color_error.size();
        ge->SetFillColor(TColor::GetColor(settings.color_error[i_spec%setting_size].c_str()));

        ge->SetFillStyle(settings.error_fill_style); // code for fill pattern

        // Add the uncertainty to the canvas if set as such
        setting_size = settings.show_error.size();
        if (settings.show_error[i_spec%setting_size])
            ge->Draw(settings.error_style.c_str());
    }

    // Draw spectra after the uncertainties so that spectrum colours are not washed out by
    // uncertainty shading. Draw in inverse order so the first plot is "at the front" of the graph.
    for (int i_spec = num_spectra-1; i_spec >= 0; i_spec--) {
        histograms[i_spec]->Draw("HIST SAME");
    }

    // Prepare the canvas for output
    TGaxis::SetMaxDigits(settings.y_digits_max); // Set maximum number of digits (i.e. for scientific notation) on yaxis
    c1->SetLogx();
    c1->Update();
    c1->Modified();
    c1->SetTickx(); // No parameter means show tick marks on both sides, labels on one
    c1->SetTicky(); // No parameter means show tick marks on both sides, labels on one
    gStyle->SetLineWidth(settings.border_width); // Set width of axis/border around the plot

    // Output the plot to file
    const char *cstr_figure_file = settings.path_output_figure.c_str();
    c1->Print(cstr_figure_file);

    return 0;
}