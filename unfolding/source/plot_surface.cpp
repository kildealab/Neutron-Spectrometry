//**************************************************************************************************
// This program reads in an arbitrary number of data series from a CSV file and plots them as a 2D
// surface plot.
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
#include "TGraph.h"
#include "TGraph2D.h"
#include "TGraphErrors.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TText.h"
#include "TVectorD.h"
#include "TVirtualPad.h"
#include "TRandom.h"

int main(int argc, char* argv[])
{
    // Read in settings
    std::string settings_file = "input/plot_surface.cfg";
    SurfaceSettings settings;
    setSurfaceSettings(settings_file, settings); // Fill settings with any user provided settings

    // Read in data
    std::vector<std::string> headers;
    std::vector<std::vector<double>> x_data;
    std::vector<std::vector<double>> y_data;
    readXYYCSV(settings.path_input_data, headers, x_data, y_data);

    int num_x = x_data[0].size();
    int num_y = y_data.size();

    // Generate the plot area
    TCanvas *c1 = new TCanvas("c1","c1",settings.x_res,settings.y_res); // Resulution of the graph (px) specified in parameters
    TGraph2D *dt = new TGraph2D();
    std::string tgraph_titles = settings.title + "; " + settings.x_label + "; " + settings.y_label + "; " + settings.z_label;
    dt->SetTitle(tgraph_titles.c_str());

    int k=0;

    // Populate with data
    for (int i=1; i<num_y; i++) {
        for (int j=1; j<num_x; j++) {
            dt->SetPoint(k,x_data[i][j],stod(headers[i]),y_data[i][j]);
            k++;
        }
    }

    // Set color scheme
    gStyle->SetPalette(settings.color_palette);
    // Set number of contours (z bins)
    gStyle->SetNumberContours(settings.num_color_bins);

    // Set plot type
    dt->Draw("COLZ");
    // dt->Draw("PCOLZ"); // 3D, no interpolation

    // Set axes ranges
    if (settings.z_min != settings.z_max) {
        dt->SetMinimum(settings.z_min);
        dt->SetMaximum(settings.z_max);
    }
    if (settings.y_min != settings.y_max){
        dt->GetYaxis()->SetLimits(settings.y_min,settings.y_max);
    }
    if (settings.x_min != settings.x_max){
        dt->GetXaxis()->SetLimits(settings.x_min,settings.x_max);
    }

    // Axes title positions
    dt->GetXaxis()->CenterTitle();
    dt->GetYaxis()->CenterTitle();
    c1->SetRightMargin(0.2); // z labels are cutoff otherwise
    dt->GetZaxis()->CenterTitle();
    dt->GetXaxis()->SetTitleOffset(1.4);
    dt->GetZaxis()->SetTitleOffset(1.8);

    // Number of ticks on axes
    if (settings.x_num_divs != 0){
        dt->GetXaxis()->SetNdivisions(settings.x_num_divs,kFALSE);
    }
    if (settings.z_num_divs != 0){
        dt->GetZaxis()->SetNdivisions(settings.z_num_divs,kFALSE);
    }

    c1->SetLogy(); // logarithmic y axis

    // Axes ticks
    c1->SetTickx(); // No parameter means show tick marks on both sides, labels on one
    c1->SetTicky(); // No parameter means show tick marks on both sides, labels on one
    gStyle->SetLineWidth(settings.border_width); // Set width of axis/border around the plot
    dt->GetZaxis()->SetTickLength(0); 

    c1->Update();
    c1->Modified();

    // Output the plot to file
    const char *cstr_figure_file = settings.path_output_figure.c_str();
    c1->Print(cstr_figure_file);
}