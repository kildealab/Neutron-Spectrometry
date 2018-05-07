#ifndef CUSTOM_CLASSES_H
#define CUSTOM_CLASSES_H

#include <stdlib.h>
#include <string>

class UnfoldingSettings {
    public:
        double norm; // vendor specfied normalization factor for the NNS used
        double error; // The target error on ratio between the experimental data points and the estimated data points from unfolding (e.g. 0.1 means the values must be within 10% of each other before the algorithm will terminate)
        double f_factor; // factor that converts measured charge (in fC) to counts per second [fA/cps]
        int cutoff; // maximum # of iterations in the unfolding algorithm
        int num_poisson_samples; // # of samples from Poisson distribution for uncertainty estimation
        // MAP specific
        double beta; // factor that affects damping of high noise spectral component in MAP method
        std::string prior;
        // Optimize specific
        int min_num_iterations;
        int max_num_iterations;
        int iteration_increment;
        double min_beta;
        double max_beta;
        std::string parameter_of_interest;

        UnfoldingSettings(); 

        void set_norm(double);
        void set_error(double);
        void set_f_factor(double);
        void set_cutoff(int);
        void set_num_poisson_samples(int);
        void set_beta(double);
        void set_prior(std::string);
        void set_min_num_iterations(int);
        void set_max_num_iterations(int);
        void set_iteration_increment(int);
        void set_min_beta(double);
        void set_max_beta(double);
        void set_parameter_of_interest(std::string);
};

#endif