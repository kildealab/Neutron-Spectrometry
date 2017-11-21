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

//==================================================================================================
// Plot a single flux spectrum (and its uncertainty) as a function of energy. Output the generated
// plot to a file using arguments passed to the function.
//==================================================================================================
int plotSpectrum(std::string figure_file_pre, std::string figure_file_suf, std::string irradiation_conditions, int num_measurements, int num_bins, std::vector<double> &energy_bins, std::vector<double> &spectrum, std::vector<double> &spectrum_uncertainty) {

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

    TCanvas *c1 = new TCanvas("c1","c1",800,600); // Resulution of the graph (px) specified in parameters
    TH1F *h1 = new TH1F("h1","h1",NBINS,edges);

    for (int i_bin = 0; i_bin < num_bins; i_bin++)
    {
        h1->Fill(energy_bins[i_bin], spectrum[i_bin]);
    }

    h1->SetStats(0);   // Do not show the stats (mean and standard deviation);
    h1->SetLineColor(kBlue);
    h1->SetLineWidth(1);
    std::ostringstream plot_title_stream;
    plot_title_stream << "Neutron flux: " << irradiation_conditions;
    std::string plot_title = plot_title_stream.str();
    h1->SetTitle(plot_title.c_str());
    h1->GetXaxis()->SetTitleOffset(1.4);
    h1->GetXaxis()->CenterTitle();
    h1->SetXTitle("Energy [MeV]");
    h1->GetYaxis()->SetTitleOffset(1.4);
    h1->GetYaxis()->CenterTitle();
    h1->SetYTitle("Fluence Rate [ncm^(-2)s^(-1)]");
    h1->Draw("HIST");  // Draw the histogram without the error bars;

    // Uncertainty plotting:

    // convert spectrum uncertainy from vector to array
    // store average bin values of adjacent bins
    // double_t s_line[NBINS]; 
    // double_t bins_line_avr[NBINS]; 
    std::vector<double> s_line;
    std::vector<double> bins_line_avr;

    for (int i_nbin = 0; i_nbin < NBINS; i_nbin++)
    {
        s_line.push_back(spectrum_uncertainty[i_nbin]);
        bins_line_avr.push_back((bins_line[i_nbin] + bins_line[i_nbin+1])/2);
    }

    // Setup plot of the uncertainty
    TGraphErrors *ge = new TGraphErrors(NBINS, &(bins_line_avr[0]), &(spectrum[0]), 0, &(s_line[0]));
    // TGraphErrors *ge = new TGraphErrors(bins_line_avr, spectrum, 0, s_line);
    ge->SetFillColor(3);
    ge->SetFillStyle(3003);
    ge->Draw("P3");

    c1->SetLogx();

    c1->Update();

    c1->Modified();

    // Output the plot to file
    std::ostringstream figure_file;
    figure_file << figure_file_pre << irradiation_conditions << figure_file_suf;
    const char *cstr_figure_file = figure_file.str().c_str();
    c1->Print(cstr_figure_file);

    return 1;
}


//==================================================================================================
// Plot a single flux spectrum (and its uncertainty) as a function of energy. Output the generated
// plot to a file using arguments passed to the function.
//==================================================================================================
int plotTwoSpectra(std::string figure_file_pre, std::string figure_file_suf, std::string irradiation_conditions, int num_measurements, int num_bins, std::vector<double> &energy_bins, std::vector<double> &spectrum, std::vector<double> &spectrum_uncertainty, std::vector<double> &spectrum2) {

    // Convert vectors to arrays for input into ROOT functions
    double ini_line[num_bins];
    double bins_line[num_bins];
    double_t edges[num_bins];
    double_t second_spectrum[num_bins];

    for (int i_bin = 0; i_bin < num_bins; i_bin++)
    {
        ini_line[i_bin] = spectrum[i_bin];
        bins_line[i_bin] = energy_bins[i_bin];
        edges[i_bin] = bins_line[i_bin];
        second_spectrum[i_bin] = spectrum2[i_bin];
    }

    // Setup plot of the spectrum
    int NBINS = num_bins-1;

    TCanvas *c1 = new TCanvas("c1","c1",800,600); // Resulution of the graph (px) specified in parameters
    TH1F *h1 = new TH1F("h1","h1",NBINS,edges);
    TH1F *h2 = new TH1F("h2","h2",NBINS,edges);

    for (int i_bin = 0; i_bin < num_bins; i_bin++)
    {
        h1->Fill(bins_line[i_bin], ini_line[i_bin]);
        h2->Fill(bins_line[i_bin], second_spectrum[i_bin]);
    }

    h2->SetStats(0);   // Do not show the stats (mean and standard deviation);
    h2->SetLineColor(kBlack);
    h2->SetLineWidth(1);

    h1->SetStats(0);   // Do not show the stats (mean and standard deviation);
    h1->SetLineColor(kRed);
    h1->SetLineWidth(1);
    std::ostringstream plot_title_stream;
    plot_title_stream << "Neutron flux: " << irradiation_conditions;
    std::string plot_title = plot_title_stream.str();
    h1->SetTitle(plot_title.c_str());
    h1->GetXaxis()->SetTitleOffset(1.4);
    h1->GetXaxis()->CenterTitle();
    h1->SetXTitle("Energy [MeV]");
    h1->GetYaxis()->SetTitleOffset(1.4);
    h1->GetYaxis()->CenterTitle();
    h1->SetYTitle("Fluence Rate [ncm^(-2)s^(-1)]");
    h2->Draw("HIST2");  // Draw the histogram without the error bars;

    h1->Draw("HIST same");  // Draw the histogram without the error bars;


    // Uncertainty plotting:

    // convert spectrum uncertainy from vector to array
    // store average bin values of adjacent bins
    double_t s_line[NBINS]; 
    double_t bins_line_avr[NBINS]; 

    for (int i_nbin = 0; i_nbin < NBINS; i_nbin++)
    {
        s_line[i_nbin] = spectrum_uncertainty[i_nbin];
        bins_line_avr[i_nbin] = (bins_line[i_nbin] + bins_line[i_nbin+1])/2;
    }

    // Setup plot of the uncertainty
    TGraphErrors *ge = new TGraphErrors(NBINS, bins_line_avr, ini_line, 0, s_line);
    ge->SetFillColor(kBlue);
    ge->SetFillStyle(3003);
    ge->Draw("P3");

    c1->SetLogx();

    c1->Update();

    c1->Modified();

    // Output the plot to file
    std::ostringstream figure_file;
    figure_file << figure_file_pre << irradiation_conditions << figure_file_suf;
    const char *cstr_figure_file = figure_file.str().c_str();
    c1->Print(cstr_figure_file);

    return 1;
}