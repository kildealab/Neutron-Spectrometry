//**************************************************************************************************
// The functions included in this module aid in performing the physics and calculations necessary
// for the neutron unfolding program.
//**************************************************************************************************

#include "physics_calculations.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <random>
#include <cmath>
#include <stdlib.h>
#include <algorithm>
#include <vector>

// mt19937 is a random number generator class based on the Mersenne Twister algorithm
// Create an mt19937 object, called mrand, that is seeded with the current time in seconds
std::mt19937 mrand(std::time(0));

//==================================================================================================
// Process Measurements: obtain sample mean measured value and standard error of sample mean for 
// each shell
//==================================================================================================
int processMeasurements(int num_measurements, int num_meas_per_shell, std::vector<double>& measurements, 
    std::vector<double>& std_errors)
{
    // Complain if mismatch with number of measurements provided and the number of measurements per
    // shell
    if (num_measurements % num_meas_per_shell != 0) {
        std::ostringstream error_message;
        error_message << "The number of measurements (" << num_measurements << ") does not evenly "
            << "divide by the number of measurements per shell (" << num_meas_per_shell << ").";
        throw std::logic_error(error_message.str());
    }

    std::vector<double> avg_measurements; // holds sample means, to replace measurements
    std::vector<double> samples; // holds the measurements comprising a given sample

    for(int i_meas=0; i_meas < num_measurements; i_meas++) {
        samples.push_back(measurements[i_meas]);

        // Once all the measurements in given sample are added to samples vector, calculate the mean
        // and standard error of the mean. Add to appropriate vectors
        if (i_meas % num_meas_per_shell == num_meas_per_shell-1) {
            double sample_mean = getMeanValueD(samples);
            avg_measurements.push_back(sample_mean);

            double std_error = getSampleMeanStandardErrorD(samples,sample_mean);
            std_errors.push_back(std_error);

            samples.clear();
        }
    }

    // Replace measurements vector with the shorter vector containing the sample means
    measurements = avg_measurements;

    return 1;
}

//==================================================================================================
// Get the mean value of a data vector containing doubles
//==================================================================================================
double getMeanValueD(std::vector<double>& data) {
    return accumulate( data.begin(), data.end(), 0.0)/data.size(); 
}


//==================================================================================================
// Get the standard error of the sample mean given a vector of data and its mean
//==================================================================================================
double getSampleMeanStandardErrorD(std::vector<double>& data, double mean) {
    int num_data = data.size();

    if (num_data == 0) {
        throw std::logic_error("Cannot calculate standard error on empty data vector.");
    }

    double numerator = 0;
    for (int i_data=0; i_data < num_data; i_data++) {
        numerator += pow((mean-data[i_data]),2);
    }
    double denominator = num_data-1;

    double single_standard_error = sqrt(numerator/denominator);

    double sample_mean_standard_error = single_standard_error/sqrt(num_data);

    return sample_mean_standard_error;
}

//==================================================================================================
// Return a 1D normalized vector of a system matrix used in MLEM-style reconstruction algorithms
//==================================================================================================
std::vector<double> normalizeResponse(int num_bins, int num_measurements, std::vector<std::vector<double>>& system_response)
{
    std::vector<double> normalized_vector;
    // Create the normalization factors to be applied to MLEM-estimated spectral values:
    //  - each element of f stores the sum of 8 elements of the transpose system (response)
    //    matrix. The 8 elements correspond to the relative contributions of each MLEM-estimated
    //    data point to the MLEM-estimated spectral value.
    for(int i_bin = 0; i_bin < num_bins; i_bin++)
    {
        double temp_value = 0;
        for(int i_meas = 0; i_meas < num_measurements; i_meas++)
        {
            temp_value += system_response[i_meas][i_bin];
        }
        normalized_vector.push_back(temp_value);
    }

    return normalized_vector;
}


//==================================================================================================
// Normalize a vector of doubles such that the largest value in the vector is set as one. The vector
// passed to this function is unaffected; a new normalized vector is returned.
//==================================================================================================
std::vector<double> normalizeVector(std::vector<double>& unnormalized_vector)
{
    std::vector<double> normalized_vector;
    double max_value = 0;
    int size = unnormalized_vector.size(); 

    // Get max value
    for (int i=0; i < size; i++) {
        if (unnormalized_vector[i] > max_value) {
            max_value = unnormalized_vector[i];
        }
    }
    // Normalize
    for (int i=0; i < size; i++) {
        normalized_vector.push_back(unnormalized_vector[i]/max_value);
    }

    return normalized_vector;
}


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
// Calculate the root-mean-square deviation of a vector of values from a "true" value.
//==================================================================================================
double calculateRMSEstimator(int size, std::vector<double> &true_vector, std::vector<double> &estimate_vector) {
    double sum_sq_diff = 0;

    // Sum the square difference of poisson sampled values from the true value
    for (int i = 0; i < size; i++) {
        sum_sq_diff += ((true_vector[i] - estimate_vector[i])*(true_vector[i] - estimate_vector[i]));
    }
    double avg_sq_diff = sum_sq_diff/size;
    double rms_diff = sqrt(avg_sq_diff);

    return rms_diff;
}


//==================================================================================================
// Calculate the normalized root mean square deviation (NRMSD) from Gaitanis et al.
//==================================================================================================
double calculateNRMSD(int size, std::vector<double> &true_vector, std::vector<double> &estimate_vector) {
    double sum_sq_diff = 0;

    // Sum the square difference of poisson sampled values from the true value
    for (int i = 0; i < size; i++) {
        sum_sq_diff += ((true_vector[i] - estimate_vector[i])*(true_vector[i] - estimate_vector[i]));
    }

    double sum_estimate = 0;
    for (int i = 0; i < size; i++) {
        sum_estimate += estimate_vector[i];
    }

    return sqrt(sum_sq_diff/sum_estimate);
}

//==================================================================================================
// Calculate chi-square from Gaitanis et al.
//==================================================================================================
double calculateChiSquaredG(int size, std::vector<double> &true_vector, std::vector<double> &estimate_vector) {
    double running_sum = 0;

    // Sum the square difference of poisson sampled values from the true value
    for (int i = 0; i < size; i++) {
        double sq_diff = ((true_vector[i] - estimate_vector[i])*(true_vector[i] - estimate_vector[i]));
        double denom = true_vector[i] + estimate_vector[i];
        running_sum += sq_diff/denom;
    }

    return 2.0*running_sum/size;
}


//==================================================================================================
// Calculate the root-mean-square deviation of a vector of values from a "true" value.
//==================================================================================================
double calculateRMSD(int num_samples, double true_value, std::vector<double> &sample_vector) {
    double sum_sq_diff = 0;

    // Sum the square difference of poisson sampled values from the true value
    for (int i_samp = 0; i_samp < num_samples; i_samp++) {
        sum_sq_diff += ((true_value - sample_vector[i_samp])*(true_value - sample_vector[i_samp]));
    }
    double avg_sq_diff = sum_sq_diff/num_samples;
    double rms_diff = sqrt(avg_sq_diff);

    return rms_diff;
}


//==================================================================================================
// Calculate the root-mean-square deviation of sampled vectors from a "true" vector.
//==================================================================================================
int calculateRMSD_vector(int num_samples, std::vector<double> &true_vector, 
    std::vector<std::vector<double>> &sampled_vectors, std::vector<double> &rms_differences)
{
    int true_vector_size = true_vector.size();

    // Loop over each element in the true vector
    for (int i = 0; i < true_vector_size; i++) {
        double sum_sq_diff = 0;
        // Sum the square difference of poisson sampled values from the true value
        for (int i_samp = 0; i_samp < num_samples; i_samp++) {
            sum_sq_diff += ((true_vector[i] - sampled_vectors[i_samp][i])*(true_vector[i] - sampled_vectors[i_samp][i]));
        }
        double avg_sq_diff = sum_sq_diff/num_samples;
        double rms_diff = sqrt(avg_sq_diff);
        rms_differences.push_back(rms_diff);
    }
    return 1;
}


//==================================================================================================
// Calculate the total ambient dose equivalent rate associated with the provided spectrum using
// weighting factors (binned by energy) provided in ICRP 74.
//==================================================================================================
double calculateDose(int num_bins, std::vector<double> &spectrum, std::vector<double> &icrp_factors) {
    double ambient_dose_eq = 0;
    int s_to_hr = 3600;
    double msv_to_psv = 1e-9;

    // Sum the dose contributions from each energy bin
    for (int i_bins=0; i_bins < num_bins; i_bins++) {
        ambient_dose_eq += spectrum[i_bins]*icrp_factors[i_bins];
    }

    // Convert from [pSv/s] to [mSv/hr]
    ambient_dose_eq *= (s_to_hr)*(msv_to_psv);

    return ambient_dose_eq;
}


//==================================================================================================
// Calculate the total measured charge associated with a series of NNS measurements
//==================================================================================================
double calculateTotalCharge(int num_measurements, std::vector<double> measurements_nc) {
    double total_charge = 0;

    for (int i_meas=0; i_meas < num_measurements; i_meas++) {
        total_charge += measurements_nc[i_meas];
    }

    return total_charge;
}


//==================================================================================================
// Calculate the total neutron flux of a neutron flux spectrum
//==================================================================================================
double calculateTotalFlux(int num_bins, std::vector<double> &spectrum) {
    double total_flux = 0;

    for (int i_bin=0; i_bin < num_bins; i_bin++) {
        total_flux += spectrum[i_bin];
    }

    return total_flux;
}


//==================================================================================================
// Calculate the total energy correction for the MAP algorithm
//==================================================================================================
double calculateTotalEnergyCorrection(std::vector<double> &energy_correction) {
    double total_energy_correction = 0;
    for (auto& n : energy_correction)
        total_energy_correction += n;

    return total_energy_correction;
}


//==================================================================================================
// Calculate the maximum deviation from 1.0 among the ratio between reconstructed CPS measurements
// measured CPS measurements
//==================================================================================================
double calculateMaxRatio(int num_measurements, std::vector<double> &mlem_ratio) {
    double max_mlem_ratio = 0;
    for (int i=0; i < num_measurements; i++) {
        if (abs(1.0-mlem_ratio[i]) > max_mlem_ratio) {
            max_mlem_ratio = abs(1.0-mlem_ratio[i]);
        }
    }
    return max_mlem_ratio;
}


//==================================================================================================
// Calculate the average deviation from 1.0 among the ratio between reconstructed CPS measurements
// measured CPS measurements
//==================================================================================================
double calculateAvgRatio(int num_measurements, std::vector<double> &mlem_ratio) {
    double avg_mlem_ratio = 0.0;
    for (int i=0; i < num_measurements; i++) {
        avg_mlem_ratio += abs(1.0-mlem_ratio[i]);
    }
    return avg_mlem_ratio/num_measurements;

}


//==================================================================================================
// Calculate the average neutron energy of a neutron flux spectrum.
// - Normalize the flux spectrum by the total flux: result is relative contribution of each energy
//      bin to the total flux
// - Multiply the energy of each bin by its relative contribution & sum
//==================================================================================================
double calculateAverageEnergy(int num_bins, std::vector<double> &spectrum, std::vector<double> &energy_bins) {
    double total_flux = calculateTotalFlux(num_bins,spectrum);
    std::vector<double> normalized_energy = spectrum;
    double avg_energy = 0;

    for (int i_bin=0; i_bin < num_bins; i_bin++) {
        avg_energy += (energy_bins[i_bin] * spectrum[i_bin]/total_flux);
    }

    return avg_energy;
}


//==================================================================================================
// Calculate the uncertainty on a sum of values.
// Implements the standard uncertainty propagation rule for sums.
//==================================================================================================
double calculateSumUncertainty(int num_values, std::vector<double> &value_uncertainties) {
    double sum_uncertainty = 0;

    for (int i=0; i < num_values; i++) {
        sum_uncertainty += pow(value_uncertainties[i],2);
    }

    return sqrt(sum_uncertainty);
}


//==================================================================================================
// Calculate the uncertainty on the average energy
// Implements standard uncertainty propagation rules for products & sums.
//==================================================================================================
double calculateEnergyUncertainty(int num_bins, std::vector<double> energy_bins, std::vector<double> spectrum, 
    std::vector<double> spectrum_uncertainty, double total_flux, double total_flux_uncertainty) 
{
    double energy_uncertainty = 0;

    for (int i_bin = 0; i_bin < num_bins; i_bin++) {
        double temp1 = energy_bins[i_bin] * spectrum[i_bin] / total_flux;
        double temp2 = sqrt(pow(spectrum_uncertainty[i_bin]/spectrum[i_bin],2)+pow(total_flux_uncertainty/total_flux,2));
        double temp3 = temp1*temp2;
        energy_uncertainty += pow(temp3,2);
    }

    energy_uncertainty = sqrt(energy_uncertainty);

    return energy_uncertainty;
}


//==================================================================================================
// Calculate the neutron source strength (shielding quanitiy of interest) of a neutron flux spectrum
// Neutron source strength = # neutrons emitted from head per Gy of photon dose delivered to 
//  isocentre

// - Empirical formula for total neutron fluence is provided in NCRP 151 pg. 42 (Eq 2.16) as a
//      function of neutron source strength. Can rearrage for source strength.
//
// Given a neutron flux spectrum:
// - Convert total flux to # neutrons per Gy per cm^2 (using dose, doserate, time, and MU->cGy calibration)
// - Divide by empirical relationship to get neutron source strength
//==================================================================================================
double calculateSourceStrength(int num_bins, std::vector<double> &spectrum, int duration, double dose_mu) {
    double total_flux = calculateTotalFlux(num_bins,spectrum);
    // Fraction of neutrons that penetrate head shielding
    double transmission_factor = 0.93; // Average of 1 for Pb and 0.85 for W
    // Surface area of treatment room [cm^2]
    double room_surface_area = 2353374.529; // For MGH.. But strength does not vary much with size
    // div by 6 = 392229 cm^2 per wall
    // l = w = 626 cm = 6.26 m (room dimension)

    // Distance from from source (Bremsstrahlung target) to point where flux was evaluated [cm]
    double distance = 100;
    // Factor to convert MU to Gy
    double mu_to_gy = 100;

    // Convert total flux to total fluence per Gy photon dose at isocenter 
    double fluence_total = total_flux * duration / dose_mu * mu_to_gy;

    // Factors in the empirical formula
    double fluence_direct_factor = transmission_factor/(4 * M_PI * pow(distance,2));
    double fluence_scatter_factor = 5.4 * transmission_factor / room_surface_area;
    double fluence_thermal_factor = 1.26 / room_surface_area;

    // Calculate source strength
    double source_strength = fluence_total / (fluence_direct_factor + fluence_scatter_factor + fluence_thermal_factor);

    return source_strength;
}

//==================================================================================================
// Calculate the indicator function, labeled as J, proposed in Bouallegue 2013 paper. This parameter
// is used in the MLEM-STOP approach as a stopping criterion when J=1.
//==================================================================================================
double calculateJFactor(int num_measurements, std::vector<double> &measurements,
    std::vector<double> &mlem_estimate) 
{
    double numerator = 0;
    double denominator = 0;
    for(int i_meas = 0; i_meas < num_measurements; i_meas++)
    {
        numerator += pow(measurements[i_meas] - mlem_estimate[i_meas],2);
        denominator += mlem_estimate[i_meas];
    }

    return numerator/denominator;
}

//==================================================================================================
// Calculate the indicator function, labeled as J2, which is a alternative metric to J. This
// indicator function also converges to 1 when within Poisson uncertainty, but is less sensitive
// to deviations from this, because is linear with respect to the MD between measurements and 
// reconstructions. 
//==================================================================================================
double calculateJFactor2(int num_measurements, std::vector<double> &measurements,
    std::vector<double> &mlem_estimate) 
{
    double numerator = 0;
    double denominator = 0;
    for(int i_meas = 0; i_meas < num_measurements; i_meas++)
    {
        numerator += abs(measurements[i_meas] - mlem_estimate[i_meas]);
        // Use of 0.798 is an emipirical finding that the expectation value of the MD of Poisson
        // samples from the mean is 0.798 * sqrt(mean) for all means (small deviations at low
        // mean values)
        denominator += 0.798*sqrt(mlem_estimate[i_meas]); 
    }

    return numerator/denominator;
}


//==================================================================================================
// Calculate the noise figure of merit in a spectrum. The noise is taken to be the relative
// standard deviation of spectral values in a region defined by the start and end energy bins
//==================================================================================================
double calculateNoise(int start_bin, int end_bin, std::vector<double>& spectrum) {
    double mean = 0.0;
    for (int i_bin=start_bin; i_bin<=end_bin; i_bin++) {
        mean += spectrum[i_bin];
    }
    mean = mean / (end_bin-start_bin+1);

    double std_dev = 0;
    for (int i_bin=start_bin; i_bin<=end_bin; i_bin++) {
        std_dev += pow((spectrum[i_bin]-mean),2);
    }
    std_dev /= (end_bin-start_bin+1);
    std_dev = sqrt(std_dev);

    return std_dev/mean;
}

//==================================================================================================
// Calculate the reduced Chi-squared parameter, promulgated by T. D. Jackson in his 2015 PhD thesis.
//==================================================================================================
double calculateChiSquared(int i_num, int num_bins, int num_measurements, std::vector<double> &spectrum, 
    std::vector<double> &measurements, std::vector<double> &mlem_ratio) 
{
    std::vector<double> mlem_estimate;

    for(int i_meas = 0; i_meas < num_measurements; i_meas++)
    {
        double temp_value = measurements[i_meas] / mlem_ratio[i_meas];
        mlem_estimate.push_back(temp_value);
    }

    double running_sum = 0;
    for(int i_meas = 0; i_meas < num_measurements; i_meas++)
    {
        double numerator = pow(measurements[i_meas] - mlem_estimate[i_meas],2);
        double denominator = pow(sqrt(measurements[i_meas]),2);
        running_sum += (numerator/denominator);
    }

    return running_sum / (num_measurements-1);
}

void calculateDerivatives(std::vector<double> &derivatives, int num_points, std::vector<int> &x_data, 
    std::vector<double> &y_data) 
{
    derivatives.clear();

    // first data point
    derivatives.push_back((y_data[2]-y_data[1])/(x_data[2]-x_data[1]));

    // Intermediate data points (not first & not last):
    for (int i=1; i < num_points-1; i++) {
        derivatives.push_back((y_data[i+1]-y_data[i-1])/(x_data[i+1]-x_data[i-1]));
    }

    // last data point
    derivatives.push_back((y_data[num_points-1]-y_data[num_points-2])/(x_data[num_points-1]-x_data[num_points-2]));
}


//==================================================================================================
// Accept a series of measurements and an estimated input spectrum and perform the MLEM algorithm
// until the true spectrum has been unfolded. Use the provided target error (error) and the maximum
// number of MLEM iterations (cutoff) to determine when to cease execution of the algorithm. Note
// that spectrum is updated as the algorithm progresses (passed by reference). Similarly for 
// mlem_ratio
//==================================================================================================
int runMLEM(int cutoff, double error, int num_measurements, int num_bins, std::vector<double> &measurements, 
    std::vector<double> &spectrum, std::vector<std::vector<double>> &nns_response, std::vector<double> &normalized_response, 
    std::vector<double> &mlem_ratio, std::vector<double> &mlem_correction, std::vector<double> &mlem_estimate) 
{
    int mlem_index; // index of MLEM iteration

    for (mlem_index = 0; mlem_index < cutoff; mlem_index++) {
        mlem_ratio.clear(); // wipe previous ratios for each iteration
        mlem_correction.clear(); // wipe previous corrections for each iteration
        mlem_estimate.clear();
        // vector that stores the MLEM-estimated data to be compared with measured data
        // std::vector<double> mlem_estimate;

        // Apply system matrix, the nns_response, to current spectral estimate to get MLEM-estimated
        // data. Save results in mlem_estimate
        // Units: mlem_estimate [cps] = nns_response [cm^2] x spectru  [cps / cm^2]
        for(int i_meas = 0; i_meas < num_measurements; i_meas++)
        {
            double temp_value = 0;
            for(int i_bin = 0; i_bin < num_bins; i_bin++)
            {
                temp_value += nns_response[i_meas][i_bin]*spectrum[i_bin];
            }
            mlem_estimate.push_back(temp_value);
        }

        // Calculate ratio between each measured data point and corresponding MLEM-estimated data point
        for(int i_meas = 0; i_meas < num_measurements; i_meas++)
        {
            mlem_ratio.push_back(measurements[i_meas]/mlem_estimate[i_meas]);
        }

        // matrix that stores the correction factor to be applied to each MLEM-estimated spectral value
        // std::vector<double> mlem_correction;

        // Create the correction factors to be applied to MLEM-estimated spectral values:
        //  - multiply transpose system matrix by ratio values
        for(int i_bin = 0; i_bin < num_bins; i_bin++)
        {
            double temp_value = 0;
            for(int i_meas = 0; i_meas < num_measurements; i_meas++)
            {
                temp_value += nns_response[i_meas][i_bin]*mlem_ratio[i_meas]/normalized_response[i_bin];
            }
            mlem_correction.push_back(temp_value);
        }

        // Apply correction factors and normalization to get new spectral estimate
        for(int i_bin=0; i_bin < num_bins; i_bin++)
        {
            spectrum[i_bin] = (spectrum[i_bin]*mlem_correction[i_bin]);
        }

        // End MLEM iterations if ratio between measured and MLEM-estimated data points is within
        // tolerace specified by 'error'
        bool continue_mlem = false;
        for (int i_meas=0; i_meas < num_measurements; i_meas++) {
            if (mlem_ratio[i_meas] >= (1+error) || mlem_ratio[i_meas] <= (1-error)) {
                continue_mlem = true;
                break;
            }   
        }
        if (!continue_mlem) {
            break;
        }
    }

    return mlem_index;
}


//==================================================================================================
// A modified version of the MLEM algorithm. A J value (Bouallegue et al 2013) is calculated at each
// iteration and unfolding is terminated when J is less than the pre-determined J threshold value.
//==================================================================================================
int runMLEMSTOP(int cutoff, int num_measurements, int num_bins, std::vector<double> &measurements, 
    std::vector<double> &spectrum, std::vector<std::vector<double>> &nns_response, std::vector<double> &normalized_response, 
    std::vector<double> &mlem_ratio, std::vector<double> &mlem_correction, std::vector<double> &mlem_estimate,
    double j_threshold, double& j_factor) 
{
    int mlem_index; // index of MLEM iteration

    for (mlem_index = 0; mlem_index < cutoff; mlem_index++) {
        mlem_ratio.clear(); // wipe previous ratios for each iteration
        mlem_correction.clear(); // wipe previous corrections for each iteration
        mlem_estimate.clear();
        // vector that stores the MLEM-estimated data to be compared with measured data
        // std::vector<double> mlem_estimate;

        // Apply system matrix, the nns_response, to current spectral estimate to get MLEM-estimated
        // data. Save results in mlem_estimate
        // Units: mlem_estimate [cps] = nns_response [cm^2] x spectru  [cps / cm^2]
        for(int i_meas = 0; i_meas < num_measurements; i_meas++)
        {
            double temp_value = 0;
            for(int i_bin = 0; i_bin < num_bins; i_bin++)
            {
                temp_value += nns_response[i_meas][i_bin]*spectrum[i_bin];
            }
            mlem_estimate.push_back(temp_value);
        }

        // Calculate ratio between each measured data point and corresponding MLEM-estimated data point
        for(int i_meas = 0; i_meas < num_measurements; i_meas++)
        {
            mlem_ratio.push_back(measurements[i_meas]/mlem_estimate[i_meas]);
        }

        // matrix that stores the correction factor to be applied to each MLEM-estimated spectral value
        // std::vector<double> mlem_correction;

        // Create the correction factors to be applied to MLEM-estimated spectral values:
        //  - multiply transpose system matrix by ratio values
        for(int i_bin = 0; i_bin < num_bins; i_bin++)
        {
            double temp_value = 0;
            for(int i_meas = 0; i_meas < num_measurements; i_meas++)
            {
                temp_value += nns_response[i_meas][i_bin]*mlem_ratio[i_meas]/normalized_response[i_bin];
            }
            mlem_correction.push_back(temp_value);
        }

        // Apply correction factors and normalization to get new spectral estimate
        for(int i_bin=0; i_bin < num_bins; i_bin++)
        {
            spectrum[i_bin] = (spectrum[i_bin]*mlem_correction[i_bin]);
        }

        // End MLEM iterations if the calculated j factor is below the threshold j value
        bool continue_mlem = true;
        j_factor = calculateJFactor(num_measurements,measurements,mlem_estimate);
        if (j_factor <= j_threshold) {
            continue_mlem = false;
        }

        if (!continue_mlem) {
            break;
        }
    }

    if (mlem_index >= cutoff && j_factor > j_threshold) {
        // Commented out 2019-11-26 b/c was printed for each poisson sample
        // std::cout << "J factor:" << j_factor << "\n";
        // std::cout << "J threshold: " << j_threshold << "\n";
        throw std::logic_error("MLEM-STOP reached cutoff # of iterations before reaching J threshold");
    }

    return mlem_index;
}


//==================================================================================================
// Calculate the J threshold for a particular set of measurements. The J threshold is taken to be
// the ratio of the average measured CPS value to the pre-determined crossover CPS value. The
// crossover value is determined graphically as explained in our manuscript.
//==================================================================================================
double determineJThreshold(int num_measurements, std::vector<double>& measurements, double cps_crossover) {
    double cps_avg = 0;
    for (int i_meas=0; i_meas < num_measurements; i_meas++) {
        cps_avg += measurements[i_meas];
    }
    cps_avg = cps_avg/num_measurements;

    double j_threshold = cps_avg/cps_crossover;
    // double j_threshold = sqrt(cps_avg/cps_crossover); // If using J2 instead of J

    return j_threshold;
}


//==================================================================================================
// Accept a series of measurements and an estimated input spectrum and perform the MLEM algorithm
// until the true spectrum has been unfolded. Use the provided target error (error) and the maximum
// number of MLEM iterations (cutoff) to determine when to cease execution of the algorithm. Note
// that spectrum is updated as the algorithm progresses (passed by reference). Similarly for 
// mlem_ratio
//==================================================================================================
int runMAP(std::vector<double> &energy_correction, double beta, std::string prior, int cutoff, double error, 
    int num_measurements, int num_bins, std::vector<double> &measurements, std::vector<double> &spectrum, 
    std::vector<std::vector<double>> &nns_response, std::vector<double> &normalized_response, std::vector<double> &mlem_ratio) 
{
    int mlem_index; // index of MLEM iteration

    for (mlem_index = 0; mlem_index < cutoff; mlem_index++) {
        mlem_ratio.clear(); // wipe previous ratios for each iteration

        // vector that stores the MLEM-estimated data to be compared with measured data
        std::vector<double> mlem_estimate;

        // Apply system matrix, the nns_response, to current spectral estimate to get MLEM-estimated
        // data. Save results in mlem_estimate
        // Units: mlem_estimate [cps] = nns_response [cm^2] x spectru  [cps / cm^2]
        for(int i_meas = 0; i_meas < num_measurements; i_meas++)
        {
            double temp_value = 0;
            for(int i_bin = 0; i_bin < num_bins; i_bin++)
            {
                temp_value += nns_response[i_meas][i_bin]*spectrum[i_bin];
            }
            mlem_estimate.push_back(temp_value);
        }

        // Calculate ratio between each measured data point and corresponding MLEM-estimated data point
        for(int i_meas = 0; i_meas < num_measurements; i_meas++)
        {
            mlem_ratio.push_back(measurements[i_meas]/mlem_estimate[i_meas]);
        }

        // matrix that stores the correction factor to be applied to each MLEM-estimated spectral value
        std::vector<double> mlem_correction;

        // Create the correction factors to be applied to MLEM-estimated spectral values:
        //  - multiply transpose system matrix by ratio values
        for(int i_bin = 0; i_bin < num_bins; i_bin++)
        {
            double temp_value = 0;
            for(int i_meas = 0; i_meas < num_measurements; i_meas++)
            {
                temp_value += nns_response[i_meas][i_bin]*mlem_ratio[i_meas];
            }
            mlem_correction.push_back(temp_value);
        }

        //------------------------------------------------------------------------------------------
        // Create the MAP energy correction factors to be incorporated in the normalization
        // std::vector<double> energy_correction;
        energy_correction.clear();

        // Quadratic prior (smoothing, no edge preservation):
        if (prior == "quadratic") {
            energy_correction.push_back(beta*pow(spectrum[0]-spectrum[1],2));
            for (int i_bin=1; i_bin < num_bins-1; i_bin++)
            {
                double temp_value = 0;
                temp_value = beta * (pow(spectrum[i_bin]-spectrum[i_bin-1],2)+pow(spectrum[i_bin]-spectrum[i_bin+1],2));          
                energy_correction.push_back(temp_value);
            }
            energy_correction.push_back(beta*pow(spectrum[num_bins-1]-spectrum[num_bins-2],2));
        }
        // Normalized quadratic prior
        else if (prior == "quadratic_normalized") {
            energy_correction.push_back(beta*sqrt(pow(spectrum[0]-spectrum[1],2))/spectrum[0]);
            for (int i_bin=1; i_bin < num_bins-1; i_bin++)
            {
                double temp_value = 0;
                temp_value = beta * 
                    sqrt(pow(spectrum[i_bin]-spectrum[i_bin-1],2)+pow(spectrum[i_bin]-spectrum[i_bin+1],2))
                    / (2*spectrum[i_bin]);          
                energy_correction.push_back(temp_value);
            }
            energy_correction.push_back(beta*sqrt(pow(spectrum[num_bins-1]-spectrum[num_bins-2],2))/spectrum[num_bins-1]);
        }
        // Median Root Prior (edge preservation by not penalizing areas of monotonic increase or decrease)
        else if (prior == "mrp") {
            int num_adjacent = 1; // on either side
            energy_correction.push_back(0); // no correction for first term
            for (int i_bin=num_adjacent; i_bin < num_bins-num_adjacent; i_bin++)
            {
                double median;

                std::vector<double> neighbours(spectrum.begin()+i_bin-num_adjacent,spectrum.begin()+i_bin+num_adjacent+1);

                // Determine if spectrum is monotonically increasing or decreasing by comparing
                // current value with its neighbours
                bool increasing = std::is_sorted(neighbours.begin(), neighbours.end());
                std::reverse(neighbours.begin(),neighbours.end());
                bool decreasing = std::is_sorted(neighbours.begin(), neighbours.end());

                // if values are monotonically increasing or decreasing, no energy correction
                if (increasing || decreasing) {
                    energy_correction.push_back(0);
                }
                // Otherwise, apply energy correction 
                else {
                    std::sort(neighbours.begin(),neighbours.end());
                    median = neighbours[num_adjacent];
                    double temp_value = beta*(spectrum[i_bin]-median)/median;
                    energy_correction.push_back(temp_value);
                }
            }
            energy_correction.push_back(0); // no correction for last term
        }
        // Custom Mean Root Prior
        else if (prior == "meanrp") {
            int num_adjacent = 1; // on either side
            energy_correction.push_back(0); // no correction for first term
            for (int i_bin=num_adjacent; i_bin < num_bins-num_adjacent; i_bin++)
            {
                double mean = 0.0;

                for (int i_n=i_bin-num_adjacent; i_n <= i_bin+num_adjacent; i_n++) {
                    mean += spectrum[i_n];
                }

                std::vector<double> neighbours(spectrum.begin()+i_bin-num_adjacent,spectrum.begin()+i_bin+num_adjacent+1);

                // Determine if spectrum is monotonically increasing or decreasing by comparing
                // current value with its neighbours
                bool increasing = std::is_sorted(neighbours.begin(), neighbours.end());
                std::reverse(neighbours.begin(),neighbours.end());
                bool decreasing = std::is_sorted(neighbours.begin(), neighbours.end());

                // if values are monotonically increasing or decreasing, no energy correction
                if (increasing || decreasing) {
                    energy_correction.push_back(0);
                }
                // Otherwise, apply energy correction 
                else {
                    double temp_value = beta*(spectrum[i_bin]-mean)/mean;
                    energy_correction.push_back(temp_value);
                }
            }
            energy_correction.push_back(0); // no correction for last term
        }
        // Custom Mean Root Prior
        else if (prior == "gaussians") {
            int num_adjacent = 1; // on either side
            energy_correction.push_back(0); // no correction for first term
            for (int i_bin=num_adjacent; i_bin < num_bins-num_adjacent; i_bin++)
            {
                double mean = 0.0;

                for (int i_n=i_bin-num_adjacent; i_n <= i_bin+num_adjacent; i_n++) {
                    mean += spectrum[i_n];
                }
                
                mean = mean / ((2*num_adjacent)+1);

                double temp_value = beta*(spectrum[i_bin]-mean)/mean;
                energy_correction.push_back(temp_value);
            }
            energy_correction.push_back(0); // no correction for last term
        }
        else {
            throw std::logic_error("Unrecognized prior: " + prior + ". Please refer to the README for allowed priors");
        }
        //------------------------------------------------------------------------------------------

        // Apply correction factors and normalization to get new spectral estimate
        for(int i_bin=0; i_bin < num_bins; i_bin++)
        {
            spectrum[i_bin] = spectrum[i_bin]*mlem_correction[i_bin]/(normalized_response[i_bin]+energy_correction[i_bin]);
        }

        // End MLEM iterations if ratio between measured and MLEM-estimated data points is within
        // tolerace specified by 'error'
        bool continue_mlem = false;
        for (int i_meas=0; i_meas < num_measurements; i_meas++) {
            if (mlem_ratio[i_meas] >= (1+error) || mlem_ratio[i_meas] <= (1-error)) {
                continue_mlem = true;
                break;
            }   
        }
        if (!continue_mlem) {
            break;
        }
    }

    return mlem_index;
}


//==================================================================================================
// Create a linearly interpolated vector of doubles with a minimum value a and a maximum value b 
// that contains N elements.
// Code retrieved from: https://gist.github.com/mortenpi/f20a93c8ed3ee7785e65
//==================================================================================================
std::vector<double> linearSpacedDoubleVector(double a, double b, std::size_t N)
{
    double h = (b - a) / static_cast<double>(N-1);
    std::vector<double> xs(N);
    std::vector<double>::iterator x;
    double val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) {
        *x = val;
    }
    return xs;
}

//==================================================================================================
// Create a linearly interpolated vector of ints with a minimum value a and a maximum value b 
// that contains N elements.
// Code retrieved from: https://gist.github.com/mortenpi/f20a93c8ed3ee7785e65
//==================================================================================================
std::vector<int> linearSpacedIntegerVector(int a, int b, std::size_t N)
{
    int h = (b - a) / (N-1);
    std::vector<int> xs(N);
    std::vector<int>::iterator x;
    int val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) {
        *x = val;
    }
    return xs;
}


//**************************************************************************************************
// RETIRED/UNUSED functions
//**************************************************************************************************
// int runMLEM_include_prev_spectrum(std::vector<double> &prev_spectrum, int cutoff, double error, 
//     int num_measurements, int num_bins, std::vector<double> &measurements, std::vector<double> &spectrum, 
//     std::vector<std::vector<double>> &nns_response, std::vector<double> &normalized_response, 
//     std::vector<double> &mlem_ratio, std::vector<double> &mlem_correction, std::vector<double> &mlem_estimate) 
// {
//     int mlem_index; // index of MLEM iteration

//     for (mlem_index = 0; mlem_index < cutoff; mlem_index++) {
//         prev_spectrum = spectrum;

//         mlem_ratio.clear(); // wipe previous ratios for each iteration
//         mlem_correction.clear(); // wipe previous corrections for each iteration
//         mlem_estimate.clear();
//         // vector that stores the MLEM-estimated data to be compared with measured data
//         // std::vector<double> mlem_estimate;

//         // Apply system matrix, the nns_response, to current spectral estimate to get MLEM-estimated
//         // data. Save results in mlem_estimate
//         // Units: mlem_estimate [cps] = nns_response [cm^2] x spectru  [cps / cm^2]
//         for(int i_meas = 0; i_meas < num_measurements; i_meas++)
//         {
//             double temp_value = 0;
//             for(int i_bin = 0; i_bin < num_bins; i_bin++)
//             {
//                 temp_value += nns_response[i_meas][i_bin]*spectrum[i_bin];
//             }
//             mlem_estimate.push_back(temp_value);
//         }

//         // Calculate ratio between each measured data point and corresponding MLEM-estimated data point
//         for(int i_meas = 0; i_meas < num_measurements; i_meas++)
//         {
//             mlem_ratio.push_back(measurements[i_meas]/mlem_estimate[i_meas]);
//         }

//         // matrix that stores the correction factor to be applied to each MLEM-estimated spectral value
//         // std::vector<double> mlem_correction;

//         // Create the correction factors to be applied to MLEM-estimated spectral values:
//         //  - multiply transpose system matrix by ratio values
//         for(int i_bin = 0; i_bin < num_bins; i_bin++)
//         {
//             double temp_value = 0;
//             for(int i_meas = 0; i_meas < num_measurements; i_meas++)
//             {
//                 temp_value += nns_response[i_meas][i_bin]*mlem_ratio[i_meas]/normalized_response[i_bin];
//             }
//             mlem_correction.push_back(temp_value);
//         }

//         // Apply correction factors and normalization to get new spectral estimate
//         for(int i_bin=0; i_bin < num_bins; i_bin++)
//         {
//             spectrum[i_bin] = (spectrum[i_bin]*mlem_correction[i_bin]);
//         }

//         // End MLEM iterations if ratio between measured and MLEM-estimated data points is within
//         // tolerace specified by 'error'
//         bool continue_mlem = false;
//         for (int i_meas=0; i_meas < num_measurements; i_meas++) {
//             if (mlem_ratio[i_meas] >= (1+error) || mlem_ratio[i_meas] <= (1-error)) {
//                 continue_mlem = true;
//                 break;
//             }   
//         }
//         if (!continue_mlem) {
//             break;
//         }
//     }

//     return mlem_index;
// }