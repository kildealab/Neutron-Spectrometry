
// This code calls the files bins.csv, respmat.csv and ini.csv and uses the MLEM algorithm to create an output matrix
// of dimensions 52*1 from a data matrix of dimensions 8*1. In addition, Poisson distribution is used to calculate a
// pseudo measured matrix of dimensions 8*1, then MLEM algorithm is used to calculate a pseudo output matrix of dimensions
// 52*1; this procedure is done 1000 times. After this, an ini_poisson matrix is created using the 1000 pseudo output matrices
// resulting in a matrix of dimensions 52*1000.
// Also, this code calculates the standard deviation and the variance by using the elements from the two matrices
// (52*1) and (52*1000).
// In addition, this code uses the ROOT libraries to plot an histogram of the neutron fluence rate (the output matrix of
// dimensions 52*1) in function of the energy bins (The result in an image.png showing the graph).

#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <ctime>
#include <cmath>

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

// mt19937 is a random number generator class based on the Mersenne Twister algorithm
// Create an mt19937 object, called mrand, that is seeded with the current time in seconds
std::mt19937 mrand(std::time(0));

// Return a poisson_distribution object, d that can be sampled from
//  - accept a parameter 'lamda', that is used as the mean & std. dev of the poisson distribution
//  - use the 'mrand' generator as the random number generator for the distribution
double poisson(double lambda)
{
    std::poisson_distribution<int> d(lambda); // initialization

    return d(mrand); // sample
}

int main()
{

    int cutoff;
    double norm;
    double error;
    int t;

    // The maximum # of iterations of the MLEM algorithm
    // Use 4000
    // Georges has compared MLEM with MAP, and found that MAP converges (reaches a plateau) at 4000,
    // whereas MLEM will continue adding noise to the signal indefinitely. 
    std::cout << "Enter the cutoff value:";
    std::cin >> cutoff;
    std::cout << "The cutoff value is equal to:" << cutoff << std::endl;

    // Normalization value for NNS measurements (vendor specified)... Need to figure out what this
    // actually is
    // Use 1.21
    std::cout << '\n';
    std::cout << "Enter the normalization value:";
    std::cin >> norm;
    std::cout << "The normalization value is equal to:" << norm << std::endl;

    // The target error on the ratio between the experimental data points and the generated data
    // points from MLEM before MLEM cutoff (e.g. 0.1 means the values must be within 10% of each
    // other before the algorithm will terminate)
    // However, the total measured fluence apparently increases with # of iterations of MLEM; which
    // renders obsolete a relative comparison in fluence between two spectra that that were generated
    // with different #s of iterations.
    // Thus use a very low target value (so the algorithm will never reach it), and will run the full
    // 4000 iterations.
    // Use 0
    std::cout << '\n';
    std::cout << "Enter the value of the error on the ratio of the experimental and processed measurements:";
    std::cin >> error;
    std::cout << "The value of the error is equal to:" << error << std::endl;

    // How long the measurement values (in nC) were acquired over (in seconds)
    // e.g. 60
    std::cout << '\n';
    std::cout << "Enter the time of measurements in seconds:";
    std::cin >> t;
    std::cout << "The value of the time is equal to:" << t << std::endl;

    double value1, value2, value3, value4, value5, value6, value7, value8;

    // Measured charge values (in nC) for each moderator & bare probe
    std::cout << '\n';
    std::cout << "Measured data (in nC):" << std::endl;
    std::cout << "Enter value1 (Bare):";
    std::cin >> value1;
    std::cout << "Enter value2 (Moderator 1):";
    std::cin >> value2;
    std::cout << "Enter value3 (Moderator 2):";
    std::cin >> value3;
    std::cout << "Enter value4 (Moderator 3):";
    std::cin >> value4;
    std::cout << "Enter value5 (Moderator 4):";
    std::cin >> value5;
    std::cout << "Enter value6 (Moderator 5):";
    std::cin >> value6;
    std::cout << "Enter value7 (Moderator 6):";
    std::cin >> value7;
    std::cout << "Enter value8 (Moderator 7):";
    std::cin >> value8;

    // measured data matrix
    // For each measured value:
    //  - multiply by normalization factor
    //  - multiply by million
    //  - divide by the measurement duration (mulitplied by 7)
    // @@Why use these numbers?
    // @@Reason for using 2D matrix?
    double data[8][1] = { {value1*norm*1000000/(7*t)} , {value2*norm*1000000/(7*t)} , {value3*norm*1000000/(7*t)} , {value4*norm*1000000/(7*t)} , {value5*norm*1000000/(7*t)} , {value6*norm*1000000/(7*t)} , {value7*norm*1000000/(7*t)} , {value8*norm*1000000/(7*t)} };


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
    // Have detector response value for each probe for each energy value (currently 52 energies)
    // Input the response from respmat.csv
    //  - values in units of [response / cm]
    //
    // The response function accounts for variable number of (n,p) reactions in He-3 for each
    // moderators, as a function of energy. Calculated by vendor using MC
    //----------------------------------------------------------------------------------------------
    int height = 8, width = 52;
    double respmat[height][width];

    std::ifstream file("respmat.csv");

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
    // Input the energies from bins.csv
    //  - values in units of [MeV]
    //  - Thi is 2D csv file, but the second value on each line is just 0 (@@not sure why this is 
    //  - necessary)
    //----------------------------------------------------------------------------------------------
    int height1 = 52, width1 = 1;
    double bins[height1][width1];

    std::ifstream file1("bins.csv");

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
    // Input from ini.csv
    //  - values are relative magnitudes (effectively fluence rates [neutrons / cm^2 s^2])
    //  - This is 2D csv file, but the second value on each line is just 0 (@@not sure why this is 
    //    necessary)
    //  - Currently (2017-08-16) input a step function (high at thermals) because:
    //    "This input spectrum was determined by selecting the spectrum that minimized the 
    //    difference between initial and reconstructed Monte Carlo neutron spectra of a Linac"
    //    (Maglieri et al. 2015)
    //    @@Why is this the case?
    //----------------------------------------------------------------------------------------------
    int height2 = 52, width2 = 1;
    double ini[height2][width2];

    std::ifstream file2("ini.csv");

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
    //std::cout << ini[height2][width2] << std::endl;

    //----------------------------------------------------------------------------------------------
    // Generate the ICRU conversion matrix (factors to convert fluence to ambient dose equivalent):
    //  - height = # of energy bins
    //  - width = @@Not sure why this is a 2D matrix
    // Input from icruconv.csv
    //  - values are @@get ICRUO 74 from Georges
    //  - This is 2D csv file, but the second value on each line is just 0 (@@not sure why this is 
    //    necessary)
    //----------------------------------------------------------------------------------------------
    int height4 = 52, width4 = 1;
    double icruconv[height4][width4];

    std::ifstream file4("icruconv.csv");

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

    for (int i = 0; i < 8; ++i)
    {
        for (int j = 0; j < 1; ++j)
        {
            std::cout << r[i][j] << ' ';
        }
        std::cout << std::endl;
    }

    //----------------------------------------------------------------------------------------------
    // Display the number of iterations of MLEM that were actually executed (<= cutoff)
    //----------------------------------------------------------------------------------------------
    std::cout << '\n';
    std::cout << "The number of iterations is equal to: " << i << std::endl;


    double data_poisson[8][1]; // analogous to 'data'. Store sampled-measured data

    int height3 = 52, width3 = 1;
    double ini_poisson[height3][width3];  // analogous to 'ini'. Store sampled-spectrum as is updated by MLEM

    double dosemat[1][1000]; // Matrix to store sampled-ambient dose equivalent values used for statistical purposes
    double poissmat[52][1000]; // Matrix to store sampled-spectra used for statistical purposes
    double dose_std = 0;

    //----------------------------------------------------------------------------------------------
    // Perform poisson analysis 1000 times
    //----------------------------------------------------------------------------------------------
    for (int m = 0; m < 1000; m++)
    {

        std::ifstream file3("ini.csv");

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
            // Converto mSv/hr
            dose_std = (3600)*(d_std)*(0.000000000001)*(1000);

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
            for (int k = 0; k < 1000; k++)
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
            sq_average[i][j] = (sq[i][j])/(1000);
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
    // multiply fluence values of the "original" MLEM-generated spectrum by ICRU conversion factors
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
    //----------------------------------------------------------------------------------------------
    double dose = (3600)*(d)*(0.000000000001)*(1000);

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
        for (int j = 0; j < 1000; j++)
            {
                sq_dose += (dosemat[i][j]-dose)*(dosemat[i][j]-dose);
            }
        }

    double sq_average_dose = sq_dose/1000;

    double s_dose = sqrt(sq_average_dose);

    std::cout << '\n';
    std::cout << "The error on the equivalent dose is: " << s_dose << " mSv/h" << std::endl;
    std::cout << '\n';


    //-------------------------------------------------------------------------------------------------------------------------------
    // ROOT plotting stuff
    //-------------------------------------------------------------------------------------------------------------------------------

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

    c1->Print("Neutron_spectrum.png");

  return 0;

}
