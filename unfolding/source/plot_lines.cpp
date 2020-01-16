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
    std::vector<std::string> headers;
    std::vector<std::vector<double>> x_data;
    std::vector<std::vector<double>> y_data;
    if (settings.data_format == "xyy") {
        readXYYCSV(settings.path_input_data, headers, x_data, y_data);
    }
    else if (settings.data_format == "xyxy") {
        readXYXYCSV(settings.path_input_data, headers, x_data, y_data);
    }
    else {
        throw std::logic_error("Unrecognized data format: " + settings.data_format);
    }

    // int num_points = x_data[0].size(); // Number of entries 
    int num_series = y_data.size();

    // Generate the plot area
    TCanvas *c1 = new TCanvas("c1","c1",settings.x_res,settings.y_res); // Resolution of the graph (px) specified in parameters
    if (settings.x_log) {
        c1->SetLogx();
    }
    if (settings.y_log) {
        c1->SetLogy();
    }
    if (settings.grayscale){
        c1->GetCanvas()->SetGrayscale();
    }

    if (settings.x_grid) {
        c1->SetGridx(1);
    }
    if (settings.y_grid) {
        c1->SetGridy(1);
    }

    // Generate the legend
    TLegend* leg = new TLegend(settings.legend_coords[0], settings.legend_coords[1], settings.legend_coords[2], 
        settings.legend_coords[3]); // with a text box
    leg->SetBorderSize(settings.legend_border_size);
    leg->SetTextSize(settings.legend_text_size);
    leg->SetFillStyle(0);

    // TVectorD xtv(num_points, &x_data[0]);
    TMultiGraph *mg = new TMultiGraph();

    // Loop through all series to be plotted
    for (int i_y = 0; i_y < num_series; i_y++) {
        // need to use TVectorD when creating TGraph:
        TVectorD xtv(x_data[i_y].size(), &x_data[i_y][0]);
        TVectorD ytv(y_data[i_y].size(), &y_data[i_y][0]); 
        TGraph* gr = new TGraph(xtv,ytv);

        // normalize all settings by the number of entries in that setting 
        // (e.g. if have 10 colors set, the 11th plot series should reuse the 1st color)
        int setting_size;

        // Styling
        setting_size = settings.plot_type.size();
        // line plots
        if (settings.plot_type[i_y%setting_size] == "l" || settings.plot_type[i_y%setting_size] == "c") {
            // color
            setting_size = settings.color_series.size();
            gr->SetLineColor(TColor::GetColor(settings.color_series[i_y%setting_size].c_str()));
            // width
            setting_size = settings.line_width.size();
            gr->SetLineWidth(settings.line_width[i_y%setting_size]);
            // style
            setting_size = settings.line_style.size();
            gr->SetLineStyle(settings.line_style[i_y%setting_size]);
        }
        // scatter plots
        else if (settings.plot_type[i_y%setting_size] == "p") {
            // // color
            setting_size = settings.color_series.size();
            gr->SetMarkerColor(TColor::GetColor(settings.color_series[i_y%setting_size].c_str()));
            // // // width
            setting_size = settings.marker_size.size();
            gr->SetMarkerSize(settings.marker_size[i_y%setting_size]);
            // // // style
            setting_size = settings.marker_style.size();
            gr->SetMarkerStyle(settings.marker_style[i_y%setting_size]);
        }
        else {
            throw std::logic_error("Unrecognized plot type: " + settings.plot_type[i_y%setting_size]);
        }

        // Add current graph to growing multigraph object
        setting_size = settings.plot_type.size();
        mg->Add(gr,settings.plot_type[i_y%setting_size].c_str());

        // Add to legend
        if (settings.legend_entries.empty()) {
            leg->AddEntry(gr, headers[i_y].c_str(), settings.plot_type[i_y%setting_size].c_str());
        }
        else {
            leg->AddEntry(gr, settings.legend_entries[i_y].c_str(), settings.plot_type[i_y%setting_size].c_str());
        }
    }

    // Add elements to canvas
    mg->Draw("a");
    if (settings.legend) {
        leg->Draw();
    }

    // Set titles
    mg->SetTitle(settings.title.c_str());
    mg->GetXaxis()->SetTitle(settings.x_label.c_str());
    mg->GetXaxis()->CenterTitle();
    mg->GetXaxis()->SetTitleOffset(settings.x_label_offset);
    mg->GetYaxis()->SetTitle(settings.y_label.c_str());
    mg->GetYaxis()->CenterTitle();
    mg->GetYaxis()->SetTitleOffset(settings.y_label_offset);

    // margins
    c1->SetLeftMargin(settings.margin_left); 
    c1->SetRightMargin(settings.margin_right); 
    c1->SetTopMargin(settings.margin_top); 
    c1->SetBottomMargin(settings.margin_bottom); 

    //Set axes properties
    if (settings.x_min != settings.x_max) {
        mg->GetXaxis()->SetLimits(settings.x_min,settings.x_max);
    }
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
        // TPaveText* pt = new TPaveText(0.15, 0.75, 0.4, 0.85, "nbNDC"); 
        // nb specifies no border, NDC specifies method of defining coordinates
        TPaveText* pt = new TPaveText(settings.textbox_coords[0], settings.textbox_coords[1], 
            settings.textbox_coords[2], settings.textbox_coords[3], "nbNDC"); 
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
    // c1->SetLogx();
    c1->Update();
    c1->Modified();
    c1->SetTickx(); // No parameter means show tick marks on both sides, labels on one
    c1->SetTicky(); // No parameter means show tick marks on both sides, labels on one
    gStyle->SetLineWidth(settings.border_width); // Set width of axis/border around the plot

    // Output the plot to file
    const char *cstr_figure_file = settings.path_output_figure.c_str();
    c1->Print(cstr_figure_file);

    return 1;
}