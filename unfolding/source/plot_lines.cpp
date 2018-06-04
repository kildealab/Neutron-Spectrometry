//**************************************************************************************************
// This program reads in an arbitrary number of data series from a CSV file and plots them as line
// graphs on a single set of axes.
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
#include "TVectorD.h"
#include "TVirtualPad.h"

int main(int argc, char* argv[])
{
    // Read in settings
    std::string settings_file = "input/plot_lines.cfg";
    PlotSettings settings;
    setPlotSettings(settings_file, settings); // Fill settings with any user provided settings

    // Read in data
    std::string input_path = settings.input_dir + settings.input_filename;
    std::vector<std::string> headers;
    std::vector<double> x_data;
    std::vector<std::vector<double>> y_data;
    readXYYCSV(input_path, headers, x_data, y_data);

    int num_points = x_data.size();
    int num_series = y_data.size();

    // Generate the plot area
    TCanvas *c1 = new TCanvas("c1","c1",settings.x_res,settings.y_res); // Resulution of the graph (px) specified in parameters

    // Generate the legend
    TLegend* leg = new TLegend(settings.legend_coords[0], settings.legend_coords[1], settings.legend_coords[2], settings.legend_coords[3]); // with a text box
    leg->SetBorderSize(0);
    leg->SetTextSize(0.035);
    leg->SetFillStyle(0);

    TVectorD xtv(num_points, &x_data[0]);
    TMultiGraph *mg = new TMultiGraph();

    // Loop through all series to be plotted
    for (int i_y = 0; i_y < num_series; i_y++) {
        TVectorD ytv(num_points, &y_data[i_y][0]); // need to use TVectorD when creating TGraph
        TGraph* gr = new TGraph(xtv,ytv);

        // Set color
        if (i_y < settings.color_series.size())        
            gr->SetLineColor(TColor::GetColor(settings.color_series[i_y].c_str()));
        // Set line width
        if (i_y < settings.line_width.size())        
            gr->SetLineWidth(settings.line_width[i_y]);
        // Set line style
        if (i_y < settings.line_style.size())        
            gr->SetLineStyle(settings.line_style[i_y]);

        // Add current graph to growing multigraph object
        if (i_y == 0)
            mg->Add(gr,"lp");
        else
            mg->Add(gr,"cp");

        // Add to legend
        if (settings.legend_entries.empty())
            leg->AddEntry(gr, headers[i_y+1].c_str(), "l"); // headers[0] is the x values header
        else
            leg->AddEntry(gr, settings.legend_entries[i_y].c_str(), "l");
    }

    // Add elements to canvas
    mg->Draw("a");
    leg->Draw();

    // Set titles
    mg->SetTitle(settings.title.c_str());
    mg->GetXaxis()->SetTitle(settings.x_label.c_str());
    mg->GetXaxis()->CenterTitle();
    mg->GetYaxis()->SetTitle(settings.y_label.c_str());
    mg->GetYaxis()->CenterTitle();

    //Set axes properties
    if (settings.x_min != settings.x_max)
        mg->GetXaxis()->SetLimits(settings.x_min,settings.x_max);
    if (settings.y_min != settings.y_max) {
        mg->SetMinimum(settings.y_min);
        mg->SetMaximum(settings.y_max);
    }
    if (settings.x_num_divs != 0)
        mg->GetXaxis()->SetNdivisions(settings.x_num_divs,kFALSE);
    if (settings.y_num_divs != 0)
        mg->GetYaxis()->SetNdivisions(settings.y_num_divs,kFALSE);

    // Optionally add a text box
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

    // Update appearance of the canvas
    c1->Update();
    c1->Modified();
    c1->SetTickx(); // No parameter means show tick marks on both sides, labels on one
    c1->SetTicky(); // No parameter means show tick marks on both sides, labels on one

    // Output the plot to file
    std::string output_path = settings.output_dir + settings.output_filename;
    const char *cstr_figure_file = output_path.c_str();
    c1->Print(cstr_figure_file);

    return 1;
}