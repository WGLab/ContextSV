#ifndef _KHMM_H
#define _KHMM_H

/// @cond
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <vector>
#include <map>
/// @endcond

// Struct for HMM (C++ RAII style)
struct CHMM
{
	int N;	// Number of states
	int M; 	// Number of observation symbols
	std::vector<std::vector<double>> A;  // Transition probability matrix
	std::vector<std::vector<double>> B;  // Emission probability matrix
	std::vector<double> pi;  // Initial state distribution
	std::vector<double> B1_mean;  // Mean of a continuous Gaussian distribution for state 1 through N
	std::vector<double> B1_sd;  // Standard deviation of B1 values, which is the same for all states
	double B1_uf;  // B1_uniform_fraction: the contribution of uniform distribution to the finite mixture model
	std::vector<double> B2_mean;  // B2_mean[1..4] is the average of B_allele_freq
	std::vector<double> B2_sd;  // B2_sd[1..4] is the standard deviation of four B_allele_freq, B2_sd[5] is specially for state1, where B is modelled as a wide normal distribution
	double B2_uf;  // B2_uniform_fraction: the fraction of uniform distribution in the finite mixture model
	int NP_flag;
	std::vector<double> B3_mean;
	std::vector<double> B3_sd;
	double B3_uf;
	int dist;
};


/************************************
*subroutines for processing continuous HMM models
************************************/

/// Read an HMM from a file
CHMM ReadCHMM (const std::string filename);

// Read a matrix
std::vector<std::vector<double>> readMatrix(std::ifstream& file, int rows, int cols);

// Read a vector
std::vector<double> readVector(std::ifstream& file, int size);

/// Run the main HMM algorithm
std::pair<std::vector<int>, double> testVit_CHMM(CHMM hmm, int T, std::vector<double>& O1, std::vector<double>& O2, std::vector<double>& pfb);

/// Viterbi algorithm
std::pair<std::vector<int>, double> ViterbiLogNP_CHMM(CHMM phmm, int T, std::vector<double>& O1, std::vector<double>& O2, std::vector<double>& pfb, double **delta, int **psi);

/// O1 emission probability (log2 coverage)
double b1iot (int state, double *mean, double *sd, double uf, double o);

/// O2 emission probability (B-allele frequency)
double b2iot (int state, double *mean, double *sd, double uf, double pfb, double b);

/// Return the probability of observing a value in a normal distribution, normalized to a range of [min_pdf, max_pdf]
double pdf_normalization(double obs, double mean, double sd);

#endif // _HMM_H_
