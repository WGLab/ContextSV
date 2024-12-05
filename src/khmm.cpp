#include "khmm.h"
#include "kc.h"

/// @cond
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
#include <stdexcept>
#include <limits>
/// @endcond

#include "utils.h"

#define STATE_CHANGE 100000.0 /*this is the expected changes (D value) in the transition matrix*/
#define VITHUGE 100000000000.0
#define FLOAT_MINIMUM 1.175494351e-38 /*this is indeed machine dependent*/
#define PROB_MAX 0.9999999999999999
#define DELTA 1

/*	This file was re-written from several subroutines from the UMDHMM package by Tapas Kanungo (Date: 15 December 1997), which has excellent framework of the implementation of Forward-Backward, Viterbi, and Baum-Welch algorithms.
	The original UMDHMM package was downloaded from http://www.kanungo.com/software/software.html. The citation for the UMDHMM program is "UMDHMM: Hidden Markov Model Toolkit," in "Extended Finite State Models of Language," A. Kornai (editor), Cambridge University Press, 1999."
	The basic framework (including variable name, subroutine name) is highly similar to the original UMDHMM package, but the actual implementation is completely different as no "discrete symbol emission" is used in PennCNV.
*/

std::pair<std::vector<int>, double> testVit_CHMM(CHMM hmm, int T, std::vector<double>& O1, std::vector<double>& O2, std::vector<double>& pfb)
{
	// T= probe, marker count
	// O1 = LRR (Log R Ratio)
	// O2 = BAF (B-Allele Freq.)
	// PFB = Genome coordinates and B-allele population frequencies
	// SNPDIST = Genome coordinates and population frequency for each SNP in the HumanHap550 array

	double **delta; // Matrix
	int **psi;		// Matrix
	delta = dmatrix(1, T, 1, hmm.N); // Allocate a TxN  double matrix (N=6 states)
	psi = imatrix(1, T, 1, hmm.N);	 // Allocate a TxN  int matrix (N=6 states)

	// Run the HMM
	std::vector<int> q;  // State sequence
	std::pair<std::vector<int>, double> state_sequence = ViterbiLogNP_CHMM(hmm, T, O1, O2, pfb, delta, psi);

	// Free the variables
	free_imatrix(psi, 1, T, 1, hmm.N);
	free_dmatrix(delta, 1, T, 1, hmm.N);

	// Pop the first element of q, which is always 0 (Done this way for 1-based
	// indexing)
	if (!state_sequence.first.empty()) {
		state_sequence.first.erase(state_sequence.first.begin());
	}

	return state_sequence;
}

// double b1iot(int state, double *mean, double *sd, double uf, double o)
double b1iot(int state, std::vector<double> mean, std::vector<double> sd, double uf, double o)
{
	// if (o < mean[1])
	// {
	// 	o = mean[1];
	// }
	// double p = uf + ((1 - uf) * pdf_normal(o, mean[state], sd[state]));

	// Get the values (0-based indexing)
	if (o < mean[0])
	{
		o = mean[0];
	}
	double p = uf + ((1 - uf) * pdf_normal(o, mean[state-1], sd[state-1]));

	return log(p);
}

// double b2iot(int state, double *mean, double *sd, double uf, double pfb, double b)
double b2iot(int state, const std::vector<double> mean, const std::vector<double> sd, double uf, double pfb, double b)
{
	// double p = 0;
	// double mean0 = mean[1];  // mean[1] = 0
	// double mean25 = mean[2];  // mean[2] = 0.25
	// double mean33 = mean[3];  // mean[3] = 0.33
	// double mean50 = mean[4];  // mean[4] = 0.5
	// double mean50_state1 = mean[5];  // mean[5] = 0.5
	// double sd0 = sd[1];  // sd[1] = 0
	// double sd25 = sd[2];  // sd[2] = 0.25
	// double sd33 = sd[3];  // sd[3] = 0.33
	// double sd50 = sd[4];  // sd[4] = 0.5
	// double sd50_state1 = sd[5];  // sd[5] = 0.5
	// p = uf;  // UF = previous alpha (transition probability)

	// Get the values (0-based indexing)
	double p = 0;
	double mean0 = mean[0];  // mean[0] = 0
	double mean25 = mean[1];  // mean[1] = 0.25
	double mean33 = mean[2];  // mean[2] = 0.33
	double mean50 = mean[3];  // mean[3] = 0.5
	double mean50_state1 = mean[4];  // mean[4] = 0.5
	double sd0 = sd[0];  // sd[0] = 0
	double sd25 = sd[1];  // sd[1] = 0.25
	double sd33 = sd[2];  // sd[2] = 0.33
	double sd50 = sd[3];  // sd[3] = 0.5
	double sd50_state1 = sd[4];  // sd[4] = 0.5
	p = uf;  // UF = previous alpha (transition probability)

	// PDF normal is the transition probability distrubution a_ij (initialized
	// as pi_n) from state i to j
	// Here, we calculate the probability of the observation symbol O_t given
	// the current state j. The observation symbol is the BAF value.
	if (state == 1)
	{
		if (b == 0)
		{
			p += (1 - uf) * cdf_normal(0, mean50_state1, sd50_state1);
		}
		else if (b == 1)
		{
			p += (1 - uf) * cdf_normal(0, mean50_state1, sd50_state1);
		}
		else
		{
			p+= (1-uf) * pdf_normal (b, mean50_state1, sd50_state1);
		}
	}
	else if (state == 2)
	{
		if (b == 0)
		{
			p += (1 - uf) * (1 - pfb) / 2;
		}
		else if (b == 1)
		{
			p += (1 - uf) * pfb / 2;
		}
		else
		{
			p+= (1-uf) * (1-pfb) * pdf_normal (b, mean0, sd0);
			p+= (1-uf) * pfb     * pdf_normal (b, 1-mean0, sd0);
		}
	}
	else if (state == 3)
	{
		if (b == 0)
		{
			p += (1 - uf) * (1 - pfb) * (1 - pfb) / 2;
		}
		else if (b == 1)
		{
			p += (1 - uf) * pfb * pfb / 2;
		}
		else
		{
			p+= (1-uf) * (1-pfb)*(1-pfb) * pdf_normal (b, mean0, sd0);
			p+= (1-uf) * 2*pfb*(1-pfb)   * pdf_normal (b, mean50, sd50);
			p+= (1-uf) * pfb*pfb         * pdf_normal (b, 1-mean0, sd0);
		}
	}
	else if (state == 4)
	{
		if (b == 0)
		{
			p += (1 - uf) * (1 - pfb) / 2;
		}
		else if (b == 1)
		{
			p += (1 - uf) * pfb / 2;
		}
		else
		{
			p+= (1-uf) * (1-pfb) * pdf_normal (b, mean0, sd0);
			p+= (1-uf) * pfb     * pdf_normal (b, 1-mean0, sd0);
		}
	}
	else if (state == 5)
	{
		if (b == 0)
		{
			p += (1 - uf) * (1 - pfb) * (1 - pfb) * (1 - pfb) / 2;
		}
		else if (b == 1)
		{
			p += (1 - uf) * pfb * pfb * pfb / 2;
		}
		else
		{
			p+= (1-uf) * (1-pfb)*(1-pfb)*(1-pfb)   * pdf_normal (b, mean0, sd0);
			p+= (1-uf) * 3*(1-pfb)*(1-pfb)*pfb     * pdf_normal (b, mean33, sd33);
			p+= (1-uf) * 3*(1-pfb)*pfb*pfb         * pdf_normal (b, 1-mean33, sd33);
			p+= (1-uf) * pfb*pfb*pfb               * pdf_normal (b, 1-mean0, sd0);
		}
	}
	else if (state == 6)
	{
		if (b == 0)
		{
			p += (1 - uf) * (1 - pfb) * (1 - pfb) * (1 - pfb) * (1 - pfb) / 2;
		}
		else if (b == 1)
		{
			p += (1 - uf) * pfb * pfb * pfb * pfb / 2;
		}
		else
		{
			p += (1-uf) * (1-pfb)*(1-pfb)*(1-pfb)*(1-pfb)   * pdf_normal (b, mean0, sd0);
			p += (1-uf) * 4*(1-pfb)*(1-pfb)*(1-pfb)*pfb     * pdf_normal (b, mean25, sd25);
			p += (1-uf) * 6*(1-pfb)*(1-pfb)*pfb*pfb         * pdf_normal (b, mean50, sd50);
			p += (1-uf) * 4*(1-pfb)*pfb*pfb*pfb             * pdf_normal (b, 1-mean25, sd25);
			p += (1-uf) * pfb*pfb*pfb*pfb                   * pdf_normal (b, 1-mean0, sd0);
		}
	}

	// Ensure that p is between FLOAT_MINIMUM and PROB_MAX
	p = std::max(FLOAT_MINIMUM, std::min(PROB_MAX, p));

	return log(p);  // Return the log probability of the observation symbol O_t
}

// PDF with normalization
double pdf_normalization(double obs, double mean, double sd)
{
	// Get the range of the PDF values
	double min_pdf = 0.0;
	double max_pdf = pdf_normal(mean, mean, sd);

	// Get the PDF value
	double pdf = pdf_normal(obs, mean, sd);

	// Normalize the PDF value
	double normalized_pdf = (pdf - min_pdf) / (max_pdf - min_pdf);
	normalized_pdf = std::max(0.0, std::min(1.0, normalized_pdf));

	return normalized_pdf;
}

std::pair<std::vector<int>, double> ViterbiLogNP_CHMM(CHMM hmm, int T, std::vector<double>& O1, std::vector<double>& O2, std::vector<double>& pfb, double **delta, int **psi)
{
	// Given the following HMM parameters:
	// - A: Transition probability matrix
	// - B: Emission probability matrix
	// - pi: Initial state distribution
	// - B1_mean: Mean of a continuous Gaussian distribution for state 1 through
	//   N
	// - B1_sd: Standard deviation of B1 values, which is the same for all
	//   states
	// - B1_uf: B1_uniform_fraction, the contribution of the uniform
	//   distribution to the finite mixture model
	// - B2_mean: Average of B_allele_freq
	// - B2_sd: Standard deviation of four B_allele_freq, B2_sd[5] is specially
	//   for state 1, where B is modelled as a wide normal distribution
	// - B2_uf: B2_uniform_fraction, the fraction of the uniform distribution in
	//   the finite mixture model
	// And the following input parameters:
	// - T: Probe, marker count (= Length of LRR array - 1)
	// - O1: LRR (Log R Ratio)
	// - O2: BAF (B-Allele Freq.)
	// - PFB: Genome coordinates and population frequency of B allele from the
	//   PFB file
	// Return the most likely state sequence Q and its likelihood using the
	// Viterbi algorithm

	int i, j; /* state indices */
	int t;	  /* time index */
	int maxvalind;
	double maxval, val;
	double **biot;
	double **A1;

	std::vector<int> q(T + 1);  // state sequence (+1 for 1-based indexing)

	// A1 is the NxN transition probability matrix from state i to j.
	// Initialize with initial values from the provided HMM file
	A1 = dmatrix(1, hmm.N, 1, hmm.N);
	for (i = 1; i <= hmm.N; i++)
	{
		for (j = 1; j <= hmm.N; j++)
		{
			// A1[i][j] = hmm.A[i][j];
			// Update for 0-based indexing
			A1[i][j] = hmm.A[i-1][j-1];
		}
	}

	/* 0. Preprocessing */

	// Pi is the initial probability distribution over states (The probability that the Markov chain will start in state i).
	// Threshold any zero values to avoid calculation issues.
	for (i = 1; i <= hmm.N; i++)
	{
		// Update to 0-based indexing
		if (hmm.pi[i-1] == 0) {
			hmm.pi[i-1] = 1e-9; /*eliminate problems with zero probability*/
		}
		hmm.pi[i-1] = log(hmm.pi[i-1]);  // Convert to log probability due to underflow
	}

	// Biot is the NxT matrix of state observation likelihoods.
	biot = dmatrix(1, hmm.N, 1, T);  // Allocate a NxT double matrix (N=6 states)
	for (i = 1; i <= hmm.N; i++)
	{
		// Loop through each probe T in the observation sequence (O1, O2), 1-based
		for (t = 1; t <= T; t++)
		{
			// Get the state observation likelihood b_j(O_t) of the observation
			// symbol O_t given the current state j

			// Calculate the O1 emission probability
			double O1_val = O1[t-1]; // Adjust for 0-based indexing

			// If there is no SNP (B-allele frequency) data, just use the LRR
			// emission probability
			if (O2[t-1] == -1)
			{

				// Calculate the O1 emission probability
				double O1_logprob = b1iot(i, hmm.B1_mean, hmm.B1_sd, hmm.B1_uf, O1_val);
				biot[i][t] = O1_logprob;

			} else {

				// Calculate the O1 emission probability
				double O1_logprob = b1iot(i, hmm.B1_mean, hmm.B1_sd, hmm.B1_uf, O1_val);

				// Calculate the O2 emission probability
				double O2_val = O2[t-1]; // Adjust for 0-based indexing
				double pfb_val = pfb[t-1]; // Adjust for 0-based indexing
				double O2_logprob = b2iot(i, hmm.B2_mean, hmm.B2_sd, hmm.B2_uf, pfb_val, O2_val);

				biot[i][t] = O1_logprob + O2_logprob;
			}
		}
	}

	/* 1. Initialization  */
	for (i = 1; i <= hmm.N; i++)
	{
		// delta[1][i] = hmm.pi[i] + biot[i][1];  // Initialize the delta matrix
		// (log probability) to the initial state distribution + the emission
		// probability
		
		// Update to 0-based indexing
		delta[1][i] = hmm.pi[i-1] + biot[i][1];  // Initialize the delta matrix
		psi[1][i] = 0;  // Initialize the psi matrix (state sequence) to 0 (no state)
	}

	/* 2. Recursion */
	// For each state j at time t, calculate the maximum probability of reaching
	// state j at time t from any of the states i at time t-1 (previous state),
	// along with observing the sequence O1, O2.
	for (t = 2; t <= T; t++)
	{
		for (j = 1; j <= hmm.N; j++)
		{
			maxval = -VITHUGE;
			maxvalind = 1;
			for (i = 1; i <= hmm.N; i++)
			{
				// Update the delta matrix (log probability) as the maximum
				// probability of being in state j at time t
				val = delta[t - 1][i] + log(A1[i][j]);
				if (val > maxval)  // Update the max value
				{
					maxval = val;
					maxvalind = i;
				}
			}


			delta[t][j] = maxval + biot[j][t];  // Update the delta matrix (log probability)
			psi[t][j] = maxvalind;  // Update the psi matrixm to store the most likely previous state at time t for state j
		}
	}

	/* 3. Termination */
	// After all observations have been processed, find the maximum probability
	// of the state sequence ending in state i at time T, along with observing
	// the sequence O1, O2.
	q[T] = 1;
	double final_lh = -VITHUGE;
	for (i = 1; i <= hmm.N; i++)
	{
		if (delta[T][i] > final_lh)
		{
			final_lh = delta[T][i];
			q[T] = i;
		}
	}

	/* 4. Path (state sequence) backtracking */
	// Backtrack through the psi matrix to find the most likely state sequence.
	// We query the psi matrix at time T using the most likely state at time T
	// from q[T]. We then backtrack through the psi matrix to find the most
	// likely state at time T-1, T-2, ..., 1.
	for (t = T - 1; t >= 1; t--)
	{
		q[t] = psi[t + 1][q[t + 1]];
	}

	// // Print t, the state, delta, biot, and psi
	// for (t = 1; t <= T; t++)
	// {
	// 	std::cout << "Time " << t << " with state " << q[t] << ":" << std::endl;
	// 	for (i = 1; i <= hmm.N; i++)
	// 	{
	// 		std::cout << "State " << i << ": delta = " << delta[t][i] << ", biot = " << biot[i][t] << ", psi = " << psi[t][i] << ", LRR = " << O1[t-1] << ", BAF = " << O2[t-1] << std::endl;
	// 	}
	// 	std::cout << std::endl;
	// }

	for (i = 1; i <= hmm.N; i++)
	{ /*recover the HMM model as original*/
		// hmm.pi[i] = exp(hmm.pi[i]);
		// Update to 0-based indexing
		hmm.pi[i-1] = exp(hmm.pi[i-1]);
	}

	free_dmatrix(biot, 1, hmm.N, 1, T);
	free_dmatrix(A1, 1, hmm.N, 1, hmm.N);

	return std::make_pair(q, final_lh);
}

CHMM ReadCHMM(const std::string filename)
{
	std::ifstream file(filename);
	if (!file.is_open())
	{
		// throw std::runtime_error("Error opening file");
		printError("Error opening file");
		return CHMM();
	}
	CHMM hmm;


	// Read M
	std::string line;
	std::getline(file, line);
	if (sscanf(line.c_str(), "M=%d", &hmm.M) != 1)
	{
		// throw std::runtime_error("Error reading M");
		printError("Error reading M");
		return CHMM();
	}

	// Read N
	std::getline(file, line);
	if (sscanf(line.c_str(), "N=%d", &hmm.N) != 1)
	{
		// throw std::runtime_error("Error reading N");
		printError("Error reading N");
		return CHMM();
	}

	// Read A
	std::getline(file, line);
	if (line != "A:")
	{
		// throw std::runtime_error("Error reading A");
		printError("Error reading A");
		return CHMM();
	}
	hmm.A = readMatrix(file, hmm.N, hmm.N);
	if (hmm.A.size() != (size_t)hmm.N || hmm.A[0].size() != (size_t)hmm.N)
	{
		// throw std::runtime_error("Error reading A");
		printError("Error reading A");
		return CHMM();
	}

	// Read B
	std::getline(file, line);
	if (line != "B:")
	{
		// throw std::runtime_error("Error reading B");
		printError("Error reading B");
		return CHMM();
	}
	hmm.B = readMatrix(file, hmm.N, hmm.M);
	if (hmm.B.size() != (size_t)hmm.N || hmm.B[0].size() != (size_t)hmm.M)
	{
		// throw std::runtime_error("Error reading B");
		printError("Error reading B");
		return CHMM();
	}

	// Read pi
	std::getline(file, line);
	if (line != "pi:")
	{
		// throw std::runtime_error("Error reading pi");
		printError("Error reading pi");
		return CHMM();
	}
	hmm.pi = readVector(file, hmm.N);
	if (hmm.pi.size() != (size_t)hmm.N)
	{
		// throw std::runtime_error("Error reading pi");
		printError("Error reading pi");
		return CHMM();
	}

	// Read B1_mean
	std::getline(file, line);
	if (line != "B1_mean:")
	{
		// throw std::runtime_error("Error reading B1_mean");
		printError("Error reading B1_mean");
		return CHMM();
	}
	hmm.B1_mean = readVector(file, hmm.N);
	if (hmm.B1_mean.size() != (size_t)hmm.N)
	{
		// throw std::runtime_error("Error reading B1_mean");
		printError("Error reading B1_mean");
		return CHMM();
	}

	// Read B1_sd
	std::getline(file, line);
	if (line != "B1_sd:")
	{
		// throw std::runtime_error("Error reading B1_sd");
		printError("Error reading B1_sd");
		return CHMM();
	}
	hmm.B1_sd = readVector(file, hmm.N);
	if (hmm.B1_sd.size() != (size_t)hmm.N)
	{
		// throw std::runtime_error("Error reading B1_sd");
		printError("Error reading B1_sd");
		return CHMM();
	}

	// Read B1_uf
	std::getline(file, line);
	if (line != "B1_uf:")
	{
		// throw std::runtime_error("Error reading B1_uf");
		printError("Error reading B1_uf");
		return CHMM();
	}
	std::getline(file, line);
	try {
		hmm.B1_uf = std::stod(line);
	} catch (const std::invalid_argument& e) {
		// throw std::runtime_error("Error reading B1_uf");
		printError("Error reading B1_uf");
		return CHMM();
	}

	// Read B2_mean
	std::getline(file, line);
	if (line != "B2_mean:")
	{
		throw std::runtime_error("Error reading B2_mean");
	}
	hmm.B2_mean = readVector(file, 5);
	if (hmm.B2_mean.size() != (size_t)5)
	{
		throw std::runtime_error("Error reading B2_mean");
	}

	// Read B2_sd
	std::getline(file, line);
	if (line != "B2_sd:")
	{
		throw std::runtime_error("Error reading B2_sd");
	}
	hmm.B2_sd = readVector(file, 5);
	if (hmm.B2_sd.size() != (size_t)5)
	{
		throw std::runtime_error("Error reading B2_sd");
	}

	// Read B2_uf
	std::getline(file, line);
	if (line != "B2_uf:")
	{
		throw std::runtime_error("Error reading B2_uf");
	}
	std::getline(file, line);
	try {
		hmm.B2_uf = std::stod(line);
	} catch (const std::invalid_argument& e) {
		throw std::runtime_error("Error reading B2_uf");
	}

	return hmm;
}

std::vector<std::vector<double>> readMatrix(std::ifstream &file, int rows, int cols)
{
    std::vector<std::vector<double>> matrix(rows, std::vector<double>(cols));
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			if (!(file >> matrix[i][j]))
			{
				throw std::runtime_error("Error reading matrix");
			}
		}
	}
	file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	return matrix;
}

std::vector<double> readVector(std::ifstream &file, int size)
{
	std::vector<double> vector(size);
	for (int i = 0; i < size; i++)
	{
		if (!(file >> vector[i]))
		{
			throw std::runtime_error("Error reading vector");
		}
	}
	file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	return vector;
}
