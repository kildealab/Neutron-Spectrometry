//**************************************************************************************************
// This code calls the files bins.csv, respmat.csv and ini.csv and uses the MLEM algorithm to create an output matrix
// of dimensions 52*1 from a data matrix of dimensions 8*1. In addition, Poisson distribution is used to calculate a
// pseudo measured matrix of dimensions 8*1, then MLEM algorithm is used to calculate a pseudo output matrix of dimensions
// 52*1; this procedure is done 1000 times. After this, an ini_poisson matrix is created using the 1000 pseudo output matrices
// resulting in a matrix of dimensions 52*1000.
// Also, this code calculates the standard deviation and the variance by using the elements from the two matrices
// (52*1) and (52*1000).
// In addition, this code uses the ROOT libraries to plot an histogram of the neutron fluence rate (the output matrix of
// dimensions 52*1) in function of the energy bins (The result in an image.png showing the graph).
//**************************************************************************************************
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <random>
#include <ctime>
#include <cmath>
#include <stdlib.h>
#include <algorithm>

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

#include <vector>

// Prototypes
double poisson(double lambda);
bool is_empty(std::ifstream& pFile);

int setfile(std::vector<std::string> &arg_vector, std::string directory, std::string arg_string, std::string default_filename, std::string &filename);
void checkUnknownParameters(std::vector<std::string> &arg_vector, std::vector<std::string> input_file_flags);
int setSettings(std::string config_file, int &cutoff, double &norm, double &error, double &f_factor, int &num_poisson_samples);
std::vector<double> getMeasurements(std::string input_file, std::string &irradiation_conditions, int &t);
int saveDose(std::string dose_file, std::string irradiation_conditions, double dose, double s_dose);
int saveSpectrum(std::string spectrum_file, std::string irradiation_conditions, double (&ini)[52][1], double (&s)[52][1], double (&bins)[52][1]);
int prepareReport(std::string report_file, std::string irradiation_conditions, std::vector<std::string> &input_files, std::vector<std::string> &input_file_flags, int cutoff, double error, double norm, double f_factor, int num_poisson_samples, std::vector<double>& data_v, int duration, std::vector<double>& bins_v, std::vector<double>& initial_v, std::vector<std::vector<double>>& response_v, int num_iterations, std::vector<double>& ratio_v, double dose, double s_dose, std::vector<double>& spectrum_v, std::vector<double>& uncertainty_v, std::vector<double>& icrp_v, std::vector<double>& subdose_v);

// Constants
std::string DOSE_HEADERS[] = {
    "Irradiation Conditions",
    "Dose Rate (mSv/hr)",
    "RMS Error (mSv/hr)"
};

std::string UNCERTAINTY_SUFFIX = "_ERROR";

// mt19937 is a random number generator class based on the Mersenne Twister algorithm
// Create an mt19937 object, called mrand, that is seeded with the current time in seconds
std::mt19937 mrand(std::time(0));

//==================================================================================================
//==================================================================================================
int main(int argc, char* argv[])
{
    // Put arguments in vector for easier processing
    std::vector<std::string> arg_vector;
    for (int i = 1; i < argc; i++) {
        arg_vector.push_back(argv[i]);
    }

    // Names of input and output directories
    std::string input_dir = "input/";
    std::string output_dir = "output/";

    // NOTE: Indices are linked between the following arrays and vectors (i.e. input_files[0]
    // corresponds to input_file_flags[0] and input_file_defaults[0])
    // Array that stores the allowed options that specify input files
    const int num_ifiles = 5;
    std::string input_file_flags_arr[num_ifiles] = {
        "--measurements",
        "--input-spectrum",
        "--energy-bins",
        "--nns-response",
        "--icrp-factors"
    };
    // Array that stores default filename for each input file
    std::string input_file_defaults_arr[num_ifiles] = {
        "measurements.txt",
        "spectrum_step.csv",
        "energy_bins.csv",
        "nns_response.csv",
        "icrp_conversions.csv"
    };

    // Convert arrays to vectors b/c easier to work with
    std::vector<std::string> input_files; // Store the actual input filenames to be used
    std::vector<std::string> input_file_flags;
    std::vector<std::string> input_file_defaults;
    for (int i=0; i<num_ifiles; i++) {
        input_files.push_back("");
        input_file_flags.push_back(input_file_flags_arr[i]);
        input_file_defaults.push_back(input_file_defaults_arr[i]);
    }

    // Use provided arguments (files) and/or defaults to determine the input files to be used
    for (int i=0; i<num_ifiles; i++) {
        setfile(arg_vector, input_dir, input_file_flags[i], input_file_defaults[i], input_files[i]);
    }

    // Notify user if unknown parameters were received
    checkUnknownParameters(arg_vector, input_file_flags);

    // Input filenames
    std::string config_file = "mlem.cfg";
    // std::string input_file;
    // std::string i_spectrum_file;
    // std::string energy_bins_file;
    // std::string response_file;
    // std::string icrp_conversions_file;

    // Output filenames
    std::string dose_file = output_dir + "output_dose.csv";
    std::string o_spectrum_file = output_dir + "output_spectra.csv"; // result (unfolded) spectrum
    std::string report_file_pre = output_dir + "REPORT_";
    std::string report_file_suf = ".txt";
    std::string figure_file_pre = output_dir + "FIGURE_";
    std::string figure_file_suf = ".png";

    // Apply some settings read in from a config file
    int cutoff; // maximum # of MLEM itereations
    double norm; // vendor specfied normalization factor for the NNS used
    double error; // The target error on ratio between the experimental data points and the estimated data points from MLEM (e.g. 0.1 means the values must be within 10% of each other before the algorithm will terminate)
    double f_factor; // factor that converts measured charge (in fC) to counts per second [fA/cps]
    int num_poisson_samples;

    bool settings_success = setSettings(config_file, cutoff, norm, error, f_factor, num_poisson_samples);
    if (!settings_success) {
        throw std::logic_error("Unable to open configuration file: " + config_file);
    }
    double f_factor_report = f_factor; // original value read in
    f_factor = f_factor / 1e6; // Convert f_factor from fA/cps to nA/cps

    // std::cout << cutoff << '\n';
    // std::cout << norm << '\n';
    // std::cout << error << '\n';
    
    // // The maximum # of iterations of the MLEM algorithm
    // // Use 4000
    // // Georges has compared MLEM with MAP, and found that MAP converges (reaches a plateau) at 4000,
    // // whereas MLEM will continue adding noise to the signal indefinitely. 
    // std::cout << "Enter the cutoff value:";
    // std::cin >> cutoff;
    // std::cout << "The cutoff value is equal to:" << cutoff << std::endl;

    // // Normalization value for NNS measurements (vendor specified)... Need to figure out what this
    // // actually is
    // // Use 1.21
    // std::cout << '\n';
    // std::cout << "Enter the normalization value:";
    // std::cin >> norm;
    // std::cout << "The normalization value is equal to:" << norm << std::endl;

    // // The target error on the ratio between the experimental data points and the generated data
    // // points from MLEM before MLEM cutoff (e.g. 0.1 means the values must be within 10% of each
    // // other before the algorithm will terminate)
    // // However, the total measured fluence apparently increases with # of iterations of MLEM; which
    // // renders obsolete a relative comparison in fluence between two spectra that that were generated
    // // with different #s of iterations.
    // // Thus use a very low target value (so the algorithm will never reach it), and will run the full
    // // 4000 iterations.
    // // Use 0
    // std::cout << '\n';
    // std::cout << "Enter the value of the error on the ratio of the experimental and processed measurements:";
    // std::cin >> error;
    // std::cout << "The value of the error is equal to:" << error << std::endl;

    // How long the measurement values (in nC) were acquired over (in seconds)
    // e.g. 60
    // std::cout << '\n';
    // std::cout << "Enter the time of measurements in seconds:";
    // std::cin >> t;
    // std::cout << "The value of the time is equal to:" << t << std::endl;

    // Read measured data (in nC) from input file
    std::string irradiation_conditions;
    int t;
    std::vector<double> data_v = getMeasurements(input_files[0], irradiation_conditions, t);

    // Convert data to existing variable structure
    // std::cout << '\n' << "------------------" << '\n';
    // convert measured charge in nC to counts per second

    double data[8][1];
    for (int index=0; index < 8; index++) {
        data[index][0] = data_v[7-index]*norm/f_factor/t;
        // std::cout <<data[index][0] << '\n';
    }

    // double value1, value2, value3, value4, value5, value6, value7, value8;

    // Measured charge values (in nC) for each moderator & bare probe
    // std::cout << '\n';
    // std::cout << "Measured data (in nC):" << std::endl;
    // std::cout << "Enter value1 (Bare):";
    // std::cin >> value1;
    // std::cout << "Enter value2 (Moderator 1):";
    // std::cin >> value2;
    // std::cout << "Enter value3 (Moderator 2):";
    // std::cin >> value3;
    // std::cout << "Enter value4 (Moderator 3):";
    // std::cin >> value4;
    // std::cout << "Enter value5 (Moderator 4):";
    // std::cin >> value5;
    // std::cout << "Enter value6 (Moderator 5):";
    // std::cin >> value6;
    // std::cout << "Enter value7 (Moderator 6):";
    // std::cin >> value7;
    // std::cout << "Enter value8 (Moderator 7):";
    // std::cin >> value8;

    // measured data matrix
    // For each measured value:
    //  - multiply by normalization factor
    //  - multiply by million
    //  - divide by the measurement duration (mulitplied by 7)
    // This converts measured values in [nC] to counts per second
    // @@Reason for using 2D matrix?
    // double data[8][1] = { {value1*norm*1000000/(7*t)} , {value2*norm*1000000/(7*t)} , {value3*norm*1000000/(7*t)} , {value4*norm*1000000/(7*t)} , {value5*norm*1000000/(7*t)} , {value6*norm*1000000/(7*t)} , {value7*norm*1000000/(7*t)} , {value8*norm*1000000/(7*t)} };


    //----------------------------------------------------------------------------------------------
    // Print out the processed measured data matrix
    //----------------------------------------------------------------------------------------------ls
    std::cout << '\n';
    std::cout << "The data matrix is equal to:" << '\n'; // newline

    //Loop over the data matrix, display each value
    for (int i = 0; i < 8; ++i)
    {
        for (int j = 0; j < 1; ++j)
        {
            std::cout << data[i][j] << ' ';
        }
        std::cout << std::endl;
    }

    //----------------------------------------------------------------------------------------------
    // Generate the detector response matrix (representing the dector response function):
    //  - height = # of measurements
    //  - width = # of energy bins
    // Detector response value for each # of moderators for each energy (currently 52 energies)
    // Input the response from input_files[3]
    //  - values in units of [cm^2]
    //
    // The response function accounts for variable number of (n,p) reactions in He-3 for each
    // moderators, as a function of energy. Calculated by vendor using MC
    //----------------------------------------------------------------------------------------------
    int height = 8, width = 52;
    double respmat[height][width];

    std::ifstream file(input_files[3]);

    //std::cout << "The response matrix is equal to:" << '\n'; // newline

    for(int row = 0; row < height; ++row)
    {
        std::string line;
        std::getline(file, line);
        if ( !file.good() )
            break;

        std::stringstream iss(line);

        for (int col = 0; col < width; ++col)
        {
            std::string val;
            std::getline(iss, val, ',');
            if ( !iss.good() )
                break;

            std::stringstream convertor(val);
            convertor >> respmat[row][col];
            //std::cout << respmat[row][col] << ' ';
        }
        //std::cout << '\n';
    }
    //std::cout << respmat[height][width] << std::endl;

    //----------------------------------------------------------------------------------------------
    // Generate the energy bins matrix:
    //  - height = # of energy bins
    //  - width = @@Not sure why this is a 2D matrix
    // Input the energies from energy bins file, input_files[2]
    //  - values in units of [MeV]
    //  - Thi is 2D csv file, but the second value on each line is just 0 (@@not sure why this is 
    //  - necessary)
    //----------------------------------------------------------------------------------------------
    int height1 = 52, width1 = 1;
    double bins[52][1];

    std::ifstream file1(input_files[2]);

    //std::cout << "The energy bins matrix is equal to:" << '\n'; // newline

    for(int row = 0; row < height1; ++row)
    {
        std::string line;
        std::getline(file1, line);
        if ( !file1.good() )
            break;

        std::stringstream iss(line);

        for (int col = 0; col < width1; ++col)
        {
            std::string val;
            std::getline(iss, val, ',');
            if ( !iss.good() )
                break;

            std::stringstream convertor(val);
            convertor >> bins[row][col];
            //std::cout << bins[row][col] << ' ';
        }
        //std::cout << '\n';
    }
    //std::cout << bins[height1][width1] << std::endl;

    //----------------------------------------------------------------------------------------------
    // Generate the inital spectrum matrix to input into MLEM algorithm:
    //  - height = # of energy bins
    //  - width = @@Not sure why this is a 2D matrix
    // Input from the input spectrum file (input_files[1])
    //  - values are neutron fluence rates [neutrons cm^-2 s^-1])
    //  - This is 2D csv file, but the second value on each line is just 0 (@@not sure why this is 
    //    necessary)
    //  - Currently (2017-08-16) input a step function (high at thermals & lower) because:
    //    "This input spectrum was determined by selecting the spectrum that minimized the 
    //    difference between initial and reconstructed Monte Carlo neutron spectra of a Linac"
    //    (Maglieri et al. 2015)
    //    @@Why is this the case?
    //----------------------------------------------------------------------------------------------
    int height2 = 52, width2 = 1;
    double ini[52][1];

    std::ifstream file2(input_files[1]);

    //std::cout << "The initial matrix is equal to:" << '\n'; // newline

    for(int row = 0; row < height2; ++row)
    {
        std::string line;
        std::getline(file2, line);
        if ( !file2.good() )
            break;

        std::stringstream iss(line);

        for (int col = 0; col < width2; ++col)
        {
            std::string val;
            std::getline(iss, val, ',');
            if ( !iss.good() )
                break;

            std::stringstream convertor(val);
            convertor >> ini[row][col];
            //std::cout << ini[row][col] << ' ';
        }
        //std::cout << '\n';
    }
    std::vector<double> initial_v;
    for (int i=0; i<52; i++) {
        initial_v.push_back(ini[i][0]);
    }
    //std::cout << ini[height2][width2] << std::endl;

    //----------------------------------------------------------------------------------------------
    // Generate the ICRP conversion matrix (factors to convert fluence to ambient dose equivalent):
    //  - height = # of energy bins
    //  - width = @@Not sure why this is a 2D matrix
    // Input from input_files[4]
    //  - values are in units of [pSv cm^2]
    //  - H values were obtained by linearly interopolating tabulated data to match energy bins used
    //  - This is 2D csv file, but the second value on each line is just 0 (@@not sure why this is 
    //    necessary)
    // Page 200 of document (ICRP 74 - ATables.pdf)
    //----------------------------------------------------------------------------------------------
    int height4 = 52, width4 = 1;
    double icruconv[height4][width4];

    std::ifstream file4(input_files[4]);

    //std::cout << "The conversion matrix is equal to:" << '\n'; // newline

    for(int row = 0; row < height4; ++row)
    {
        std::string line;
        std::getline(file4, line);
        if ( !file4.good() )
            break;

        std::stringstream iss(line);

        for (int col = 0; col < width4; ++col)
        {
            std::string val;
            std::getline(iss, val, ',');
            if ( !iss.good() )
                break;

            std::stringstream convertor(val);
            convertor >> icruconv[row][col];
            //std::cout << icruconv[row][col] << ' ';
        }
        //std::cout << '\n';
    }
    //std::cout << icruconv[height4][width4] << std::endl;

    //----------------------------------------------------------------------------------------------
    // Run the MLEM algorithm, iterating <cutoff> times.
    // Final result, i.e. MLEM-estimated spectrum, outputted in 'ini' matrix
    //----------------------------------------------------------------------------------------------
    int i; // index of MLEM iteration
    int j; 
    int k;

    // matrix that stores the ratio between measured data points and corresponding MLEM data point
    double r[8][1]; 

    for (i = 1.0; i < cutoff; i++) {

        // matrix that stores the MLEM-estimated data to be compared with measured data
        double d[8][1];

        // Apply system matrix 'respmat' to current spectral estimate 'ini' to get MLEM-estimated
        // data. Save results in 'd'
        // Units: d [cps] = respmat [cm^2] x ini [cps / cm^2]
        for(int i0 = 0; i0 < 8; i0++)
        {
            for(j = 0; j < 1; j++)
            {
            d[i0][j]=0;
                for(k = 0; k < 52; k++)
                {
                d[i0][j]=d[i0][j]+(respmat[i0][k]*ini[k][j]);
                }
            }
        }

        // Create the transpose matrix of the system matrix (for upcoming backprojection)
        double trans_respmat[52][8];

        for(int i1 = 0; i1 < 8; ++i1)
        for(j = 0; j < 52; ++j)
        {
           trans_respmat[j][i1]=respmat[i1][j];
        }

        // std::cout << "The transpose matrix of respmat is equal to:" << '\n'; // newline

        // for (int i2 = 0; i2 < 52; ++i2)
        // {
        //     for (int j = 0; j < 8; ++j)
        //     {
        //         std::cout << trans_respmat[i2][j] << ' ';
        //     }
        //     std::cout << std::endl;
        // }

        // Calculate ratio between each measured data point and corresponding MLEM-estimated data point
        for(int i3 = 0; i3 < 8; i3++)
        {
            for(j = 0; j < 1; j++)
            {
            r[i3][j]=0;
            r[i3][j]=r[i3][j]+(data[i3][j]/d[i3][j]);
            }
        }

        // matrix that stores the correction factor to be applied to each MLEM-estimated spectral value
        double c[52][1];

        // Create the correction factors to be applied to MLEM-estimated spectral values:
        //  - multiply transpose system matrix by ratio values
        for(int i4 = 0; i4 < 52; i4++)
        {
            for(j = 0; j < 1; j++)
            {
            c[i4][j]=0;
                for(k = 0; k < 8; k++)
                {
                c[i4][j]=c[i4][j]+(trans_respmat[i4][k]*r[k][j]);
                }
            }
        }

        double one[8][1] = { {1} , {1} , {1} , {1} , {1} , {1} , {1} , {1} };

        // matrix that stores the normalization factors (sensitivity) to be applied to each MLEM-estimated
        // spectral value
        double f[52][1];

        // Create the normalization factors to be applied to MLEM-estimated spectral values:
        //  - each element of f stores the sum of 8 elements of the transpose system (response)
        //    matrix. The 8 elements correspond to the relative contributions of each MLEM-estimated
        //    data point to the MLEM-estimated spectral value.
        for(int i5 = 0; i5 < 52; i5++)
        {
            for(j = 0; j < 1; j++)
            {
            f[i5][j]=0;
                for(k = 0; k < 8; k++)
                {
                f[i5][j]=f[i5][j]+(trans_respmat[i5][k]*one[k][j]);
                }
            }
        }


        //Temporary matrix to store @@unnormalized corrected estimated spectral values
        double ini1[52][1];
        // Apply correction factors to previous estimate of spectral values
        for(int i6 = 0; i6 < 52; i6++)
        {
            for(j = 0; j < 1; j++)
            {
            ini1[i6][j]=0;
            ini1[i6][j]=ini1[i6][j]+(ini[i6][j]*c[i6][j]);
            }
        }

        // Apply normalization factors to corrected estimate of spectral values
        for(int i7 = 0; i7 < 52; i7++)
        {
            for(j = 0; j < 1; j++)
            {
            ini[i7][j]=0;
            ini[i7][j]=ini[i7][j]+(ini1[i7][j]/f[i7][j]);
            }
        }

        // Prematurely end MLEM iterations if ratio between measured and MLEM-estimated data points
        // is within tolerace specified by 'error'
        if ( r[0][0] < (1+error) && r[0][0] > (1-error) && r[1][0] < (1+error) && r[1][0] > (1-error) && r[2][0] < (1+error) && r[2][0] > (1-error) && r[3][0] < (1+error) && r[3][0] > (1-error) && r[4][0] < (1+error) && r[4][0] > (1-error) && r[5][0] < (1+error) && r[5][0] > (1-error) && r[6][0] < (1+error) && r[6][0] > (1-error) && r[7][0] < (1+error) && r[7][0] > (1-error) )
        {
            break;
        }

    }

    //----------------------------------------------------------------------------------------------
    // Display the result (output) matrix of the MLEM algorithm, which represents reconstructed 
    // spectral data.
    //----------------------------------------------------------------------------------------------
    std::cout << '\n';
    std::cout << "The output matrix is equal to:" << '\n'; // newline

    for (int i8 = 0; i8 < 52; ++i8)
    {
        for (int j = 0; j < 1; ++j)
        {
            std::cout << ini[i8][j] << ' ';
        }
        std::cout << std::endl;
    }

    //----------------------------------------------------------------------------------------------
    // Display the ratio (error) matrix of the MLEM algorithm, which represents the deviation of the
    // MLEM-generated measured charge values (which correspond to the above MLEM-generated measured 
    // spectrum) from the actual measured charge values. 
    //----------------------------------------------------------------------------------------------
    std::cout << '\n';
    std::cout << "The error matrix is equal to:" << '\n'; // newline
    // double avg_ratio = 0; // Calculate average error in ratio across 8 measuremnts
    // double max_ratio = 0; // maximum error in ratio across 8 measurements

    for (int i = 0; i < 8; ++i)
    {
        for (int j = 0; j < 1; ++j)
        {
            std::cout << r[i][j] << ' ';
            // double myval = std::abs(r[i][j]-1.0);
            // avg_ratio += myval;
            // if (myval > max_ratio) {
            //     max_ratio = myval;
            // }
        }
        std::cout << std::endl;
    }

    //----------------------------------------------------------------------------------------------
    // Display the number of iterations of MLEM that were actually executed (<= cutoff)
    //----------------------------------------------------------------------------------------------
    std::cout << '\n';
    std::cout << "The number of iterations is equal to: " << i << std::endl;
    int num_interations = i;


    double data_poisson[8][1]; // analogous to 'data'. Store sampled-measured data

    int height3 = 52, width3 = 1;
    double ini_poisson[height3][width3];  // analogous to 'ini'. Store sampled-spectrum as is updated by MLEM

    double dosemat[1][num_poisson_samples]; // Matrix to store sampled-ambient dose equivalent values used for statistical purposes
    double poissmat[52][num_poisson_samples]; // Matrix to store sampled-spectra used for statistical purposes
    double dose_std = 0;

    //----------------------------------------------------------------------------------------------
    // Perform poisson analysis 1000 times
    //----------------------------------------------------------------------------------------------
    for (int m = 0; m < num_poisson_samples; m++)
    {

        std::ifstream file3(input_files[1]);

        //std::cout << "The initial poisson matrix is equal to:" << '\n'; // newline

        //------------------------------------------------------------------------------------------
        // Generate the inital spectrum matrix to input into poisson MLEM algorithm:
        //  - height = # of energy bins
        //  - width = @@Not sure why this is a 2D matrix
        // @@Seems like this is the same matrix as imported above into 'ini' (should be able to just
        // make a copy before running MLEM above, rather than re-import from csv here)
        //------------------------------------------------------------------------------------------
        for(int row = 0; row < height3; ++row)
        {
            std::string line;
            std::getline(file3, line);
            if ( !file3.good() )
                break;

            std::stringstream iss(line);

            for (int col = 0; col < width3; ++col)
            {
                std::string val;
                std::getline(iss, val, ',');
                if ( !iss.good() )
                    break;

                std::stringstream convertor(val);
                convertor >> ini_poisson[row][col];
                //std::cout << ini[row][col] << ' ';
            }
            //std::cout << '\n';
        }
        //std::cout << ini_poisson[height3][width3] << std::endl;

        //------------------------------------------------------------------------------------------
        // sample from poission distribution created for each measured data point, to create a new
        // set of sampled-data points. Store in 'data_poisson'
        //------------------------------------------------------------------------------------------
        for (int j = 0; j < 1; j++)
        {
            for (int i = 0; i <8; i++)
            {
                data_poisson[i][j] = poisson(data[i][j]);
            }
        }

        //std::cout << '\n';
        //std::cout << "The " << (m+1) << "th poisson data matrix is equal to:" << '\n'; // newline

        //for (int i1 = 0; i1 < 8; ++i1)
        //{
            //for (int j1 = 0; j1 < 1; ++j1)
            //{
                //std::cout << data_poisson[i1][j1] << ' ';
            //}
            //std::cout << std::endl;
        //}

        //------------------------------------------------------------------------------------------
        // Do MLEM as before
        //------------------------------------------------------------------------------------------
        for (int i = 1; i < cutoff; i++) {

            double d_poisson[8][1]; // MLEM-estimated data

            // Apply system matrix to estimated spectrum
            for(int i0 = 0; i0 < 8; i0++)
            {
                for(j = 0; j < 1; j++)
                {
                d_poisson[i0][j]=0;
                    for(k = 0; k < 52; k++)
                    {
                        d_poisson[i0][j]=d_poisson[i0][j]+(respmat[i0][k]*ini_poisson[k][j]);
                    }
                }
            }

            // Create transpose system matrix
            double trans_respmat[52][8];
            for(int i1 = 0; i1 < 8; ++i1)
            for(j = 0; j < 52; ++j)
            {
               trans_respmat[j][i1]=respmat[i1][j];
            }

            //std::cout << "The transpose matrix of respmat is equal to:" << '\n'; // newline

            //for (int i2 = 0; i2 < 52; ++i2)
            //{
                //for (int j = 0; j < 8; ++j)
                //{
                    //std::cout << trans_respmat[i2][j] << ' ';
                //}
                //std::cout << std::endl;
            //}

            double r_poisson[8][1]; // ratio between measured and MLEM-estimated data

            // calculate ratios
            for(int i3 = 0; i3 < 8; i3++)
            {
                for(j = 0; j < 1; j++)
                {
                    r_poisson[i3][j]=0;
                    r_poisson[i3][j]=r_poisson[i3][j]+(data_poisson[i3][j]/d_poisson[i3][j]);
                }
            }

            double c_poisson[52][1]; // correction factors

            // calculate correction factors
            for(int i4 = 0; i4 < 52; i4++)
            {
                for(j = 0; j < 1; j++)
                {
                    c_poisson[i4][j]=0;
                    for(k = 0; k < 8; k++)
                    {
                        c_poisson[i4][j]=c_poisson[i4][j]+(trans_respmat[i4][k]*r_poisson[k][j]);
                    }
                }
            }

            double one[8][1] = { {1} , {1} , {1} , {1} , {1} , {1} , {1} , {1} };

            double f_poisson[52][1]; // normalization factors

            // calculate normalization factors
            for(int i5 = 0; i5 < 52; i5++)
            {
                for(j = 0; j < 1; j++)
                {
                f_poisson[i5][j]=0;
                    for(k = 0; k < 8; k++)
                    {
                        f_poisson[i5][j]=f_poisson[i5][j]+(trans_respmat[i5][k]*one[k][j]);
                    }
                }
            }

            double ini1_poisson[52][1]; // temporary spectral matrix

            // apply correction factors to previous spectral estimate
            for(int i6 = 0; i6 < 52; i6++)
            {
                for(j = 0; j < 1; j++)
                {
                    ini1_poisson[i6][j]=0;
                    ini1_poisson[i6][j]=ini1_poisson[i6][j]+(ini_poisson[i6][j]*c_poisson[i6][j]);
                }
            }

            // apply normalization factors to corrected spectral estimate
            for(int i7 = 0; i7 < 52; i7++)
            {
                for(j = 0; j < 1; j++)
                {
                    ini_poisson[i7][j]=0;
                    ini_poisson[i7][j]=ini_poisson[i7][j]+(ini1_poisson[i7][j]/f_poisson[i7][j]);
                }
            }

            // multiply fluence values of the estimated spectrum by ICRU conversion factors
            // Converts from [n cm^-2 s^-1] to [pSv s^-1]
            double multiplication_std[52][1];
            for (int i = 0; i < 52; i++)
            {
                for (int j = 0; j < 1; j++)
                {
                    multiplication_std[i][j] = ini_poisson[i][j]*icruconv[i][j];
                }
            }

            // Get the total dose, as contributed to by fluence value in each energy bin
            double d_std = 0;
            for (int i = 0; i < 52; i++)
            {

                for (int j = 0; j < 1; j++)
                {
                    d_std += multiplication_std[i][j];
                }
            }
            // Convert from [pSv/s] to [mSv/hr]
            dose_std = (d_std)*(3600)*(1e-9);

            // Prematurely end MLEM iterations if ratio between measured and MLEM-estimated data points
            // is within tolerace specified by 'error'
            if ( r_poisson[0][0] < (1+error) && r_poisson[0][0] > (1-error) && r_poisson[1][0] < (1+error) && r_poisson[1][0] > (1-error) && r_poisson[2][0] < (1+error) && r_poisson[2][0] > (1-error) && r_poisson[3][0] < (1+error) && r_poisson[3][0] > (1-error) && r_poisson[4][0] < (1+error) && r_poisson[4][0] > (1-error) && r_poisson[5][0] < (1+error) && r_poisson[5][0] > (1-error) && r_poisson[6][0] < (1+error) && r_poisson[6][0] > (1-error) && r_poisson[7][0] < (1+error) && r_poisson[7][0] > (1-error) )
            {
                break;
            }
        } // Finish MLEM

        //std::cout << '\n';
        //std::cout << "The " << (m+1) << "th output poisson matrix is equal to:" << '\n'; // newline

        //for (int i8 = 0; i8 < 52; ++i8)
        //{
            //for (int j = 0; j < 1; ++j)
            //{
                //std::cout << ini_poisson[i8][j] << ' ';
            //}
            //std::cout << std::endl;
        //}

        // Transfer the MLEM-estimated spectrum into the growing 'poissmat' (which stores the
        // spectrum generated for each statistical iteration [not each iteration of MLEM])
        for (int j = 0; j < 1; j++)
        {
            for (int i = 0; i <52; i++)
            {
                poissmat[i][m] = ini_poisson[i][j];
            }
        }

        // Transfer the calculated dose value into the growing 'dosemat' (whic stores the dose value
        // calculated for each statistical iteration)
        for (int i = 0; i < 1; i++)
        {
            dosemat[i][m] = dose_std;
        }

    } // Finish Poisson looping

    //std::cout << '\n';
    //std::cout << "The dose matrix is equal to:" << '\n'; // newline

    //for (int i = 0; i < 1; ++i)
    //{
        //for (int j = 0; j < 1000; ++j)
        //{
            //std::cout << dosemat[i][j] << ' ';
        //}
        //std::cout << std::endl;
    //}

    //----------------------------------------------------------------------------------------------
    // Calculate the summed squared difference of all the sampled spectra from the "original" MLEM-
    // generated spectrum
    //----------------------------------------------------------------------------------------------
    double sq[52][1]; 
    for (int i = 0; i < 52; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            for (int k = 0; k < num_poisson_samples; k++)
            {
                sq[i][j] += (poissmat[i][k]-ini[i][j])*(poissmat[i][k]-ini[i][j]);
            }
        }
    }

    //----------------------------------------------------------------------------------------------
    // Calculate the average squared difference of a sampled spectra from the "original" MLEM-
    // generated spectrum
    //----------------------------------------------------------------------------------------------
    double sq_average[52][1];
    for (int i =0; i < 52; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            sq_average[i][j] = (sq[i][j])/(num_poisson_samples);
        }
    }

    //----------------------------------------------------------------------------------------------
    // Calculate the RMS difference of a sampled spectra from the "original" MLEM-generated spectrum
    //----------------------------------------------------------------------------------------------
    double s[52][1];
    for (int i =0; i < 52; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            s[i][j] = sqrt(sq_average[i][j]);
        }
    }

    //----------------------------------------------------------------------------------------------
    // Print the RMS difference matrix
    //
    // @@misleading output here; is not the standard deviation matrix. Is the RMS difference matrix
    //----------------------------------------------------------------------------------------------
    std::cout << '\n'; // newline
    std::cout << "The standard deviation matrix is equal to:" << '\n'; // newline

    for (int i = 0; i < 52; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            std::cout << s[i][j] << ' ';
        }
        std::cout << std::endl;
    }

    //----------------------------------------------------------------------------------------------
    // multiply fluence values of the "original" MLEM-generated spectrum by ICRP conversion factors
    //----------------------------------------------------------------------------------------------
    double multiplication[52][1];

    for (int i = 0; i < 52; i++)
    {
        for (int j = 0; j < 1; j++)
        {
            multiplication[i][j] = ini[i][j]*icruconv[i][j];
        }
    }

    //----------------------------------------------------------------------------------------------
    // Get the total dose, as contributed to by fluence value in each energy bin
    //----------------------------------------------------------------------------------------------
    double d = 0;

    for (int i = 0; i < 52; i++)
    {

        for (int j = 0; j < 1; j++)
        {
            d += multiplication[i][j];
        }
    }

    //----------------------------------------------------------------------------------------------
    // Converto mSv/hr
    // Display result
    // Convert from [pSv/s] to [mSv/hr]
    //----------------------------------------------------------------------------------------------
    double dose = (d)*(3600)*(1e-9);

    std::cout << '\n';
    std::cout << "The equivalent dose is: " << dose << " mSv/h" << std::endl;
    std::cout << '\n';

    //----------------------------------------------------------------------------------------------
    // Get RMS difference of sampled dose from MLEM-estimated dose
    // Display result
    //----------------------------------------------------------------------------------------------
    double sq_dose = 0;

    for (int i = 0; i < 1; i++)
        {
        for (int j = 0; j < num_poisson_samples; j++)
            {
                sq_dose += (dosemat[i][j]-dose)*(dosemat[i][j]-dose);
            }
        }

    double sq_average_dose = sq_dose/num_poisson_samples;

    double s_dose = sqrt(sq_average_dose);

    std::cout << '\n';
    std::cout << "The error on the equivalent dose is: " << s_dose << " mSv/h" << std::endl;
    std::cout << '\n';

    //----------------------------------------------------------------------------------------------
    // Save results to file
    //----------------------------------------------------------------------------------------------
    saveDose(dose_file, irradiation_conditions, dose, s_dose);
    std::cout << "Saved calculated dose to " << dose_file << "\n";
    saveSpectrum(o_spectrum_file, irradiation_conditions, ini, s, bins);
    std::cout << "Saved unfolded spectrum to " << o_spectrum_file << "\n";
    //************************************
    // Convert arrays to vectors
    std::vector<double> bins_v;
    std::vector<std::vector<double>> response_v;
    std::vector<double> ratio_v;
    std::vector<double> spectrum_v;
    std::vector<double> uncertainty_v;
    std::vector<double> icrp_v;
    std::vector<double> subdose_v;
    for (int i=0; i < 52; i++) {
        bins_v.push_back(bins[i][0]);
        std::vector<double> new_column;
        for (int j=0; j<8; j++) {
            new_column.push_back(respmat[j][i]);
        }
        response_v.push_back(new_column);
        spectrum_v.push_back(ini[i][0]);
        uncertainty_v.push_back(s[i][0]);
        icrp_v.push_back(icruconv[i][0]);
        subdose_v.push_back(multiplication[i][0]*(3600)*(1e-9));
    }
    for (int i=0; i<8; i++) {
        ratio_v.push_back(r[i][0]);
    }
    //************************************
    std::string report_file = report_file_pre + irradiation_conditions + report_file_suf;
    prepareReport(report_file, irradiation_conditions, input_files, input_file_flags, cutoff, error, norm, f_factor_report, num_poisson_samples, data_v, t, bins_v, initial_v, response_v, num_interations, ratio_v, dose, s_dose, spectrum_v, uncertainty_v, icrp_v, subdose_v);
    std::cout << "Generated summary report: " << report_file << "\n\n";
    // std::cout << "Avg ratio: " << avg_ratio/8.0 << "\n";
    // std::cout << "Max ratio: " << max_ratio << "\n\n";

    //----------------------------------------------------------------------------------------------
    // ROOT plotting stuff
    //----------------------------------------------------------------------------------------------

    //-----------------------------------------------------------------------------------------------
    // Creating the line matrices for (ini and bins matrices) to be used in the plotting
    //-----------------------------------------------------------------------------------------------

    double ini_line[52];

    for (int i = 0; i < 52; i++)
    {
        for (int j = 0; j < 1; j++)
        {
        ini_line[i] = ini[i][j];
        }
    }

    //std::cout << '\n';
    //std::cout << "The output line matrix is equal to:" << '\n'; // newline

    //for (int i = 0; i < 52; ++i)
    //{
    //std::cout << ini_line[i] << ' ';
    //}

    double bins_line[52];

    for (int i = 0; i < 52; i++)
    {
        for (int j = 0; j < 1; j++)
        {
        bins_line[i] = bins[i][j];
        }
    }

    //std::cout << '\n';
    //std::cout << "The energy bins line matrix is equal to:" << '\n'; // newline

    //for (int i = 0; i < 52; ++i)
    //{
    //std::cout << bins_line[i] << ' ';
    //}


    //---------------------------------------------------------------------------------------------
    // ROOT plotting procedure
    //---------------------------------------------------------------------------------------------

    int NBINS = 51;

    double_t edges[NBINS + 1];

    for (int i = 0; i < 52; i++)
    {
        edges[i] = bins_line[i];
    }

    TCanvas *c1 = new TCanvas("c1","c1",800,600); // Resulution of the graph is 800*600 pixels.

    TH1F *h1 = new TH1F("h1","h1",NBINS,edges);

    for (int i = 0; i < 52; i++)
    {
        h1->Fill(bins_line[i], ini_line[i]);
    }

    h1->SetStats(0);   // Do not show the stats (mean and standard deviation);
    h1->SetLineColor(kBlue);
    h1->SetLineWidth(1);
    h1->SetTitle("NEUTRON SPECTRUM");
    h1->GetXaxis()->SetTitleOffset(1.4);
    h1->GetXaxis()->CenterTitle();
    h1->SetXTitle("Energy [MeV]");
    h1->GetYaxis()->SetTitleOffset(1.4);
    h1->GetYaxis()->CenterTitle();
    h1->SetYTitle("Fluence Rate [ncm^(-2)s^(-1)]");
    h1->Draw("HIST");  // Draw the histogram without the error bars;


    double_t s_line[51];

    for (int i = 0; i < 51; i++)
    {
        for (int j = 0; j < 1; j++)
        {
        s_line[i] = s[i][j];
        }
    }

    //std::cout << '\n';
    //std::cout << "The standard deviation line matrix is equal to:" << '\n'; // newline

    //for (int i = 0; i < 51; ++i)
    //{
    //std::cout << s_line[i] << ' ';
    //}

    double_t bins_line_avr[51];

    for (int i = 0; i < 51; i++)
    {
    bins_line_avr[i] = (bins_line[i] + bins_line[i+1])/2;
    }

    //std::cout << '\n';
    //std::cout << "The average energy bins line matrix is equal to:" << '\n'; // newline

    //for (int i = 0; i < 51; ++i)
    //{
    //std::cout << bins_line_avr[i] << ' ';
    //}

    double_t ini_line_e[51];

    for (int i = 0; i < 51; i++)
    {
        for (int j = 0; j < 1; j++)
        {
        ini_line_e[i] = ini[i][j];
        }
    }

    //std::cout << '\n';
    //std::cout << "The output line matrix is equal to:" << '\n'; // newline

    //for (int i = 0; i < 51; ++i)
    //{
    //std::cout << ini_line_e[i] << ' ';
    //}

    TGraphErrors *ge = new TGraphErrors(51, bins_line_avr, ini_line_e, 0, s_line);
    ge->SetFillColor(3);
    ge->SetFillStyle(3003);
    ge->Draw("P3");

    c1->SetLogx();

    c1->Update();

    c1->Modified();

    std::ostringstream figure_file;
    figure_file << figure_file_pre << irradiation_conditions << figure_file_suf;
    const char *cstr_figure_file = figure_file.str().c_str();
    c1->Print(cstr_figure_file);

  return 0;

}

//**************************************************************************************************
// Helper functions
//**************************************************************************************************

//==================================================================================================
// Return a poisson_distribution object, d that can be sampled from
//  - accept a parameter 'lamda', that is used as the mean & std. dev of the poisson distribution
//  - use the 'mrand' generator as the random number generator for the distribution
//==================================================================================================
double poisson(double lambda)
{
    std::poisson_distribution<int> d(lambda); // initialization

    return d(mrand); // sample
}

//==================================================================================================
// Determine if a file is empty
//==================================================================================================
bool is_empty(std::ifstream& pFile)
{
    return pFile.peek() == std::ifstream::traits_type::eof();
}

//==================================================================================================
// Check arguments passed to the main function for a matching flag, and assign value for the input
// filename
// Args:
//  - arg_vector: Vector containing all arguments passed to the main function (at command line)
//  - directory: Directory (path) at which the input file is located (e.g. input/)
//  - arg_string: The argument string that indicates what the file is (e.g. --measurements)
//  - default_file: The default filename to be used for the arg_string provided
//  - filename: The variable (passed by reference) that will be assigned a value, representing the
//      the filename to be used
//==================================================================================================
int setfile(std::vector<std::string> &arg_vector, std::string directory, std::string arg_string, std::string default_filename, std::string &filename) {
    // Place an iterator at an index of the vector, where the matching arg_string was found
    std::vector<std::string>::iterator iter_args = std::find(arg_vector.begin(), arg_vector.end(), arg_string);

    // If a matching argument was found, ensure that a subsequent argument (i.e. filename to use) 
    // was provided. Apply it if so. If not, throw an error
    if( iter_args != arg_vector.end()) {
        iter_args += 1;
        if (iter_args != arg_vector.end()) {
            filename = directory + *iter_args;
        }
        else {
            // throw error saying no file provided for *iter_measurements -= 1
            iter_args -= 1;
            std::cout << "Error: no file provided for argument: " << *iter_args << "\n";
        }
    }
    // If no match was found for the target arg_string within arg_vector, use the default filename
    // i.e. The user did not specify which file to use
    else {
        filename = directory + default_filename;
    }

    return 1;
}

//==================================================================================================
// Identify any options that were provided to the main function that are not supported and notify
// the user.
//==================================================================================================
void checkUnknownParameters(std::vector<std::string> &arg_vector, std::vector<std::string> input_file_flags) {
    for (int i=0; i<arg_vector.size(); i++) {
        std::vector<std::string>::iterator iter_args = std::find(input_file_flags.begin(), input_file_flags.end(), arg_vector[i]);
        if(iter_args == input_file_flags.end()) {
            std::cout << "Warning: Ignored unknown argument " << arg_vector[i] << "\n";
        }
        else {
            i += 1;
        }
    }
   
}

//==================================================================================================
// Retrieve settings from a configuration file 'config_file' and save values in relevant variables
//==================================================================================================
int setSettings(std::string config_file, int &cutoff, double &norm, double &error, double &f_factor, int &num_poisson_samples) {
    std::ifstream cfile(config_file);
    std::string line;
    std::vector<double> settings;

    // If file is able to be read
    if (cfile.is_open())
    {
        // loop through each line in the file, extract the value for each setting into 'token'
        while ( getline (cfile,line) )
        {
            // settings format: setting_name=value
            std::string delimiter = "=";
            std::string token = line.substr(line.find(delimiter)+1); // substring from '=' to end of string
            // std::cout << token << '\n';
            settings.push_back(atof(token.c_str())); // convert str to double, insert into vector
        }
        cfile.close();
    }
    // Problem opening the file
    else {
        return false;
    }

    // Assign the settings values to function parameters passed by reference
    cutoff = (int)settings[0];
    norm = settings[1];
    error = settings[2];
    f_factor = settings[3];
    num_poisson_samples = settings[4];
    return true;
}

//==================================================================================================
// Read measurement data from file
// data_vector is ordered from index:0 storing the value for 7 moderators to index:7 storing the value
// for 0 moderators
//==================================================================================================
std::vector<double> getMeasurements(std::string input_file, std::string &irradiation_conditions, int &t) {
    std::ifstream ifile(input_file);
    if (!ifile.is_open()) {
        //throw error
        std::cout << "Unable to open input file:" + input_file + '\n';
    }

    // Load header information from 'ifile'
    getline(ifile,irradiation_conditions);
    // removes: carriage return '\r' from the string (which causes weird string overwriting)
    irradiation_conditions.erase( std::remove(irradiation_conditions.begin(), irradiation_conditions.end(), '\r'), irradiation_conditions.end() );

    // Extract measurement duration
    std::string t_string;
    getline(ifile,t_string);
    t = atoi(t_string.c_str());

    // Loop through file, get measurement data
    std::string line;
    std::vector<double> data_vector;
    while (getline(ifile,line)) {
        std::istringstream line_stream(line);
        std::string stoken; // store individual values between delimiters on a line

        // Loop through each line, delimiting at commas
        while (getline(line_stream, stoken, ',')) {
            data_vector.push_back(atof(stoken.c_str())); // add data to the vector
        }
    }

    // print elements of data vector
    // int vector_size = data_vector.size();
    // for (int i=0; i<vector_size; i++) {
    //     std::cout << data_vector[i] << '\n';
    // }
    ifile.close();
    std::cout << "Data successfully retrieved from " + input_file + '\n';
    return data_vector;
}

//==================================================================================================
// Save calculated dose (and its error) to file
//==================================================================================================
int saveDose(std::string dose_file, std::string irradiation_conditions, double dose, double s_dose) {
    // determine if file exists
    std::ifstream checkfile(dose_file);
    bool file_empty = is_empty(checkfile);
    checkfile.close();

    // If file does not exist, create it. If file exists, open it
    std::ofstream dfile;
    dfile.open(dose_file, std::ios_base::app);

    // Add header line to start of file, if file was empty
    if (file_empty) {
        std::ostringstream header_stream;
        header_stream << DOSE_HEADERS[0] << "," << DOSE_HEADERS[1] << "," << DOSE_HEADERS[2] << "\n";
        std::string header = header_stream.str();
        dfile << header; 
    }

    // Append new line of data
    std::ostringstream new_data_stream;
    new_data_stream << irradiation_conditions << "," << dose << "," << s_dose << "\n";
    std::string new_data = new_data_stream.str();
    dfile << new_data; 

    dfile.close();
    return 1;
}

//==================================================================================================
// Save calculated spectrum (and error spectrum) to file
//==================================================================================================
int saveSpectrum(std::string spectrum_file, std::string irradiation_conditions, double (&ini)[52][1], double (&s)[52][1], double (&bins)[52][1]) {
    // determine if file exists
    std::ifstream sfile(spectrum_file);
    bool file_empty = is_empty(sfile);
    bool file_exists = sfile.good();

    std::vector<std::string> sfile_lines;

    // If the file exists and is not empty: append 2 new columns to existing rows (strings)
    // First new column: the spectrum, Second new column: the uncertainty on the spectrum
    if (file_exists && !file_empty) {
        // Retrieve and update the header
        int index = 0;
        std::string header;
        getline(sfile, header);
        header = header + "," + irradiation_conditions + "," + irradiation_conditions + UNCERTAINTY_SUFFIX + "\n";
        sfile_lines.push_back(header);

        // Loop through existing lines in file (including header)
        // For each, string-concatenate the new data
        // Store new lines in vector 'sfile_lines'
        std::string line;
        while (getline(sfile,line)) {
            //Note: using stream b/c have to convert doubles to strings
            std::ostringstream line_stream;
            line_stream << line << "," << ini[index][0] << "," << s[index][0] << "\n";
            line = line_stream.str();
            sfile_lines.push_back(line);
            index++;
        }
    }
    // If the file does not exist or was empty: create each line for the file
    // Column 1: energy bins, Column 2: the spectrum, Column 3: the uncertainty on spectrum
    else {
        // Prepare header
        std::string header;
        header = "Energy (MeV)," + irradiation_conditions + "," + irradiation_conditions + UNCERTAINTY_SUFFIX + "\n";
        sfile_lines.push_back(header);

        // Prepare contents
        std::string line;
        int spectrum_length = sizeof ini / sizeof ini[0];
        for (int index=0; index<spectrum_length; index++) {
            //Note: using stream b/c have to convert doubles to strings
            std::ostringstream line_stream;
            line_stream << bins[index][0] << "," << ini[index][0] << "," << s[index][0] << "\n";
            line = line_stream.str();
            sfile_lines.push_back(line);
        }
    }
    sfile.close();

    // Rewrite file using lines stored in vector 'sfile_lines'
    std::ofstream nfile(spectrum_file);
    int vector_size = sfile_lines.size();
    for (int i=0; i<vector_size; i++) {
        nfile << sfile_lines[i];
    }

    return 1;
}

//==================================================================================================
// Generate a textfile report of pertinent information from execution of this program. The contents
// of this function are separated by headers indicating the type of information printed to the
// report in the the corresponding section.
//==================================================================================================
int prepareReport(std::string report_file, std::string irradiation_conditions, std::vector<std::string> &input_files, std::vector<std::string> &input_file_flags, int cutoff, double error, double norm, double f_factor, int num_poisson_samples, std::vector<double>& data_v, int duration, std::vector<double>& bins_v, std::vector<double>& initial_v, std::vector<std::vector<double>>& response_v, int num_iterations, std::vector<double>& ratio_v, double dose, double s_dose, std::vector<double>& spectrum_v, std::vector<double>& uncertainty_v, std::vector<double>& icrp_v, std::vector<double>& subdose_v) {
    std::string HEADER_DIVIDE = "************************************************************************************************************************\n";
    std::string SECTION_DIVIDE = "\n========================================================================================================================\n\n";
    std::string COLSTRING = "--------------------";
    int sw = 30; // settings column width
    int cw = 20; // data column width
    int rw = 9; // NNS response column width

    std::ofstream rfile(report_file);

    //----------------------------------------------------------------------------------------------
    // Header
    //----------------------------------------------------------------------------------------------
    rfile << HEADER_DIVIDE;
    rfile << "Neutron Spectrometry Report\n\n";
    rfile << std::left << std::setw(sw) << "Irradiation Specs: " << irradiation_conditions << "\n";
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    rfile << std::left << std::setw(sw) << "Date report was generated: " << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << "\n";
    rfile << "Input arguments (files) used:\n";
    for (int i=0; i<input_files.size(); i++) {
        std::string tempstring = "    " + input_file_flags[i];
        rfile << std::left <<std::setw(sw) << tempstring << input_files[i] << "\n";
    }
    rfile << HEADER_DIVIDE << "\n";

    //----------------------------------------------------------------------------------------------
    // Settings
    //----------------------------------------------------------------------------------------------
    rfile << "Settings\n\n";
    rfile << std::left << std::setw(sw) << "MLEM max # of iterations:" << cutoff << "\n";
    rfile << std::left << std::setw(sw) << "MLEM target ratio:" << error << "\n";
    rfile << std::left << std::setw(sw) << "NNS normalization factor:" << norm << "\n";
    rfile << std::left << std::setw(sw) << "NNS calibration factor:" << f_factor << " fA/cps\n";
    rfile << std::left << std::setw(sw) << "Number of poisson samples:" << num_poisson_samples << "\n";
    rfile << SECTION_DIVIDE;

    //----------------------------------------------------------------------------------------------
    // Measurement
    //----------------------------------------------------------------------------------------------
    rfile << "Measured Data (measurement duration: " << duration << "s)\n\n";
    rfile << std::left << std::setw(cw) << "# of moderators" << "Charge (nC)\n";
    rfile << std::left << std::setw(cw) << COLSTRING << COLSTRING << "\n";
    for (int i=0; i<data_v.size(); i++) {
        rfile << std::left << std::setw(cw) << i << data_v[i] << "\n";
    }
    rfile << SECTION_DIVIDE;

    //----------------------------------------------------------------------------------------------
    // Inputs
    //----------------------------------------------------------------------------------------------
    rfile << "Inputs (Number of energy bins: " << bins_v.size() << ")\n\n";
    rfile << std::left << std::setw(cw) << "Energy bins" << std::setw(cw) << "Input spectrum" << "| NNS Response by # of moderators (cm^2)\n";
    rfile << std::left << std::setw(cw) << "(MeV)" << std::setw(cw) << "(n cm^-2 s^-1)" << "| ";
    for (int j=0; j<response_v[0].size(); j++) {
        rfile << std::left << std::setw(rw) << j;
    }
    rfile << "\n";
    rfile << std::left << std::setw(cw) << COLSTRING << std::setw(cw) << COLSTRING << "--";
    for (int j=0; j<response_v[0].size(); j++) {
        rfile << "---------";
    }
    rfile << "\n";

    for (int i=0; i<bins_v.size(); i++) {
        rfile << std::left << std::setw(cw) << bins_v[i] << std::setw(cw) << initial_v[i] << "| ";
        for (int j=0; j<response_v[0].size(); j++) {
            rfile << std::left << std::setw(rw) << response_v[i][j];
        }
        rfile << "\n";
    }
    rfile << SECTION_DIVIDE;

    //----------------------------------------------------------------------------------------------
    // MLEM Processing
    //----------------------------------------------------------------------------------------------
    rfile << "MLEM information\n\n";
    rfile << std::left << std::setw(sw) << "# of iterations: " << num_iterations << "/" << cutoff << "\n\n";
    rfile << "Final MLEM ratio = measured charge / estimated charge:\n";
    int thw = 13; // NNS response column width
    //row 1
    rfile << std::left << std::setw(thw) << "# moderators" << "| ";
    for (int j=0; j<ratio_v.size(); j++) {
        rfile << std::left << std::setw(rw) << j;
    }
    rfile << "\n";
    // row 2
    rfile << std::left << std::setw(thw) << "-------------|-";
    for (int j=0; j<response_v[0].size(); j++) {
        rfile << "---------";
    }
    rfile << "\n";
    // row thw
    rfile << std::left << std::setw(thw) << "ratio" << "| ";
    for (int j=0; j<ratio_v.size(); j++) {
        rfile << std::left << std::setw(rw) << ratio_v[j];
    }
    rfile << "\n";
    rfile << SECTION_DIVIDE;

    //----------------------------------------------------------------------------------------------
    // Results
    //----------------------------------------------------------------------------------------------
    rfile << "Results\n\n";
    rfile << std::left << std::setw(sw) << "Ambient dose equivalent:" << dose << " mSv/hr\n";
    rfile << std::left << std::setw(sw) << "Uncertainty:" << s_dose << " mSv\n\n";
    rfile << std::left << std::setw(cw) << "Energy bins" << std::setw(cw) << "Unfolded spectrum" << std::setw(cw) << "Uncertainty" << std::setw(cw) << "| ICRP H factor" << "Ambient Dose Equiv.\n";
    rfile << std::left << std::setw(cw) << "(MeV)" << std::setw(cw) << "(n cm^-2 s^-1)" << std::setw(cw) << "(n cm^-2 s^-1)" << std::setw(cw) << "| (pSv/cm^2)" << "(mSv/hr)\n";;
    rfile << std::left << std::setw(cw) << COLSTRING << std::setw(cw) << COLSTRING << std::setw(cw) << COLSTRING << std::setw(cw) << COLSTRING << COLSTRING << "\n";
    for (int i=0; i<bins_v.size(); i++) {
        std::ostringstream icrp_string;
        icrp_string << "| " <<icrp_v[i];
        rfile << std::left << std::setw(cw) << bins_v[i] << std::setw(cw) << spectrum_v[i] << std::setw(cw) << uncertainty_v[i] << std::setw(26) << icrp_string.str() << subdose_v[i] << "\n";
    }

    rfile.close();
    return 1;
}
