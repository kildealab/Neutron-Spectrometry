//**************************************************************************************************
// The functions included in this module aid in producing plots using the ROOT library.
//**************************************************************************************************

#include "root_helpers.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <stdlib.h>
#include <vector>

// Root
#include "TAxis.h"
#include "TGaxis.h"
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
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TText.h"

//==================================================================================================
// Plot a single flux spectrum (and its uncertainty) as a function of energy. Output the generated
// plot to a file using arguments passed to the function.
//==================================================================================================
int plotSpectrum(std::string path_figure, std::string irradiation_conditions, 
    int num_measurements, int num_bins, std::vector<double> &energy_bins, std::vector<double> &spectrum, 
    std::vector<double> &spectrum_uncertainty_upper,std::vector<double> &spectrum_uncertainty_lower) 
{

    // Convert vectors to arrays for input into ROOT functions
    double ini_line[num_bins];
    double bins_line[num_bins];
    double_t edges[num_bins];

    for (int i_bin = 0; i_bin < num_bins; i_bin++)
    {
        ini_line[i_bin] = spectrum[i_bin];
        bins_line[i_bin] = energy_bins[i_bin];
        edges[i_bin] = bins_line[i_bin];
    }

    // Setup plot of the spectrum
    int NBINS = num_bins-1;

    TCanvas *c1 = new TCanvas("c1","c1",2400,1800); // Resulution of the graph (px) specified in parameters
    TH1F *h1 = new TH1F("h1","h1",NBINS,edges);

    for (int i_bin = 0; i_bin < num_bins; i_bin++)
    {
        h1->Fill(energy_bins[i_bin], spectrum[i_bin]);
    }

    h1->SetStats(0);   // Do not show the stats (mean and standard deviation);
    h1->SetLineColor(kBlack);
    h1->SetLineWidth(10);
    std::ostringstream plot_title_stream;
    plot_title_stream << "Neutron fluence spectrum: " << irradiation_conditions;
    std::string plot_title = plot_title_stream.str();
    h1->SetTitle(plot_title.c_str());
    h1->GetXaxis()->SetTitleOffset(1.4);
    h1->GetXaxis()->CenterTitle();
    h1->SetXTitle("Energy (MeV)");
    h1->GetYaxis()->SetTitleOffset(1.4);
    h1->GetYaxis()->CenterTitle();
    h1->SetYTitle("Fluence (n #upoint cm^{-2} s^{-1})");
    h1->SetTickLength(0.015,"xy"); // Length of tick marks (default = 0.02)
    h1->Draw("HIST");  // Draw the histogram without the error bars;

    // Uncertainty plotting:

    // convert spectrum uncertainy from vector to array
    // store average bin values of adjacent bins
    // double_t s_line[NBINS]; 
    // double_t bins_line_avr[NBINS]; 
    std::vector<double> s_line;
    std::vector<double> bins_line_avr;

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
        // error_upper_array.push_back(spectrum_uncertainty_upper[i_bin]);
        // error_lower_array.push_back(spectrum_uncertainty_lower[i_bin]);
    }

    // for (int i_nbin = 0; i_nbin < NBINS; i_nbin++)
    // {
    //     s_line.push_back(spectrum_uncertainty[i_nbin]);
    //     bins_line_avr.push_back((bins_line[i_nbin] + bins_line[i_nbin+1])/2);
    // }

    // Setup plot of the uncertainty
    // TGraphErrors *ge = new TGraphErrors(NBINS, &(bins_line_avr[0]), &(spectrum[0]), 0, &(s_line[0]));

    TGraphAsymmErrors *ge = new TGraphAsymmErrors(
        num_bins, &(bins_log_avg[0]), &(spectrum[0]), &(xerror_lower_array[0]), &(xerror_upper_array[0]), 
        &(spectrum_uncertainty_lower[0]), &(spectrum_uncertainty_upper[0])
    );
    // TGraphErrors *ge = new TGraphErrors(bins_line_avr, spectrum, 0, s_line);
    ge->SetFillColor(920);
    ge->SetFillStyle(3001);
    ge->Draw("P2");
    h1->Draw("HIST SAME");  // Draw the histogram without the error bars;

    TGaxis::SetMaxDigits(3); // Set maximum number of digits (i.e. for scientific notation) on yaxis

    c1->SetLogx();

    c1->Update();

    c1->Modified();

    c1->SetTickx(); // No parameter means show tick marks on both sides, labels on one
    c1->SetTicky(); // No parameter means show tick marks on both sides, labels on one

    gStyle->SetLineWidth(5);

    // Output the plot to file
    const char *cstr_figure_file = path_figure.c_str();
    c1->Print(cstr_figure_file);
    return 0;
}