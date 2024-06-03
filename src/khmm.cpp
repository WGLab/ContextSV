#include "khmm.h"
#include "kc.h"

/// @cond
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <iomanip>
/// @endcond

#define STATE_CHANGE 100000.0 /*this is the expected changes (D value) in the transition matrix*/
#define VITHUGE 100000000000.0
#define FLOAT_MINIMUM 1.175494351e-38 /*this is indeed machine dependent*/
#define PROB_MAX 0.9999999999999999
#define DELTA 1

/*	This file was re-written from several subroutines from the UMDHMM package by Tapas Kanungo (Date: 15 December 1997), which has excellent framework of the implementation of Forward-Backward, Viterbi, and Baum-Welch algorithms.
	The original UMDHMM package was downloaded from http://www.kanungo.com/software/software.html. The citation for the UMDHMM program is "UMDHMM: Hidden Markov Model Toolkit," in "Extended Finite State Models of Language," A. Kornai (editor), Cambridge University Press, 1999."
	The basic framework (including variable name, subroutine name) is highly similar to the original UMDHMM package, but the actual implementation is completely different as no "discrete symbol emission" is used in PennCNV.
*/

// Entry point
// std::vector<int> testVit_CHMM(CHMM hmm, int T, double *O1, double *O2, double *pfb, int *snpdist, double *plogproba)
std::pair<std::vector<int>, double> testVit_CHMM(CHMM hmm, int T, std::vector<double>& O1, std::vector<double>& O2, std::vector<double>& pfb)
{
	// T= probe, marker count (= Length of LRR array - 1)
	// O1 = LRR (Log R Ratio)
	// O2 = BAF (B-Allele Freq.)
	// PFB = Genome coordinates and population frequency of B allele from the
	// PFB file
	// Note: PFB is replaced with Q (state sequence) after running the HMM
	// SNPDIST = Genome coordinates and population frequency for each SNP in the HumanHap550 array
	// PLOGPROBA = LOGPROB, log probability for each of the 6 states
	// Note: PLOGPROBA is replaced with DELTA (delta matrix) after running the
	// HMM, which is used to calculate the probability of the most likely state

	double **delta; // Matrix
	int **psi;		// Matrix
	delta = dmatrix(1, T, 1, hmm.N); // Allocate a TxN  double matrix (N=6 states)
	psi = imatrix(1, T, 1, hmm.N);	 // Allocate a TxN  int matrix (N=6 states)
	
	// Set initial posterior log probabilities for each state
	std::vector<double>	plogproba(hmm.N + 1, -VITHUGE);

	// Run the HMM
	std::vector<int> q;  // State sequence
	std::pair<std::vector<int>, double> state_sequence = ViterbiLogNP_CHMM(hmm, T, O1, O2, pfb, delta, psi, plogproba);
	// q = ViterbiLogNP_CHMM(hmm, T, O1, O2, pfb, delta, psi, plogproba);
	// q = state_sequence.first;

	// Use the forward algorithm to calculate the probability of the most likely
	// state sequence


	// Free the variables
	free_imatrix(psi, 1, T, 1, hmm.N);
	free_dmatrix(delta, 1, T, 1, hmm.N);

	// Pop the first element of q, which is always 0 (Done this way for 1-based
	// indexing)
	state_sequence.first.erase(state_sequence.first.begin());
	// q.erase(q.begin());

	return state_sequence;
}

// Emission probability: Calculate the state observationn likelihood b_j(O_t) of
// the observation symbol O_t given the current state j
// O_t is the LRR value
double b1iot(int state, double *mean, double *sd, double uf, double o)
{
	// double p = uf;  // UF = Uniform distribution term

	// Get the normalized PDF value
	double normalized_pdf = pdf_normalization(o, mean[state], sd[state]);

	// // Calculate the PDF term for the observation symbol O_t given the
	// // current state j
	// double pdf_value = pdf_normal(o, mean[state], sd[state]);

	// // Normalize the PDF value
	// double normalized_pdf = (pdf_value - min_pdf) / (max_pdf - min_pdf);
	// normalized_pdf = std::max(0.0, std::min(1.0, normalized_pdf));

	// Update the emission probability
	// p += (1 - uf) * normalized_pdf;
	double p = uf + ((1 - uf) * normalized_pdf);

	// PDF normal is the transition probability distrubution a_ij (initialized
	// as pi_n) from state i to j
	// TODO: Sometimes the PDF normal is greater than 1, which is incorrect
	// p += (1 - uf) * pdf_normal(o, mean[state], sd[state]);	

	// // Throw an error if p is not between 0 and 1
	// if (p < 0 || p > 1)
	// {
	// 	std::cerr << "ERROR [3] = " << p << std::endl;
	// 	std::cerr << "Normalized O = " << normalized_pdf << std::endl;
	// 	std::cerr << "UF = " << uf << std::endl;
	// }

	// Prevent divide by zero error
	// if (p == 0)
	// {
	// 	p = FLOAT_MINIMUM;
	// }

	// Ensure that p is between FLOAT_MINIMUM and PROB_MAX
	p = std::max(FLOAT_MINIMUM, std::min(PROB_MAX, p));

	// Return the log probability
	return log(p);
}

// Emission probability: Calculate the state observationn likelihood b_j(O_t) of
// the observation symbol O_t given the current state j
// O_t is the BAF value
double b2iot(int state, double *mean, double *sd, double uf, double pfb, double b)
{
	double p = 0;
	double mean0 = mean[1];  // mean[1] = 0
	double mean25 = mean[2];  // mean[2] = 0.25
	double mean33 = mean[3];  // mean[3] = 0.33
	double mean50 = mean[4];  // mean[4] = 0.5
	double mean50_state1 = mean[5];  // mean[5] = 0.5
	double sd0 = sd[1];  // sd[1] = 0
	double sd25 = sd[2];  // sd[2] = 0.25
	double sd33 = sd[3];  // sd[3] = 0.33
	double sd50 = sd[4];  // sd[4] = 0.5
	double sd50_state1 = sd[5];  // sd[5] = 0.5

	p = uf;  // UF = previous alpha (transition probability)

	// PDF normal is the transition probability distrubution a_ij (initialized
	// as pi_n) from state i to j
	// Here, we calculate the probability of the observation symbol O_t given
	// the current state j. The observation symbol is the BAF value.
	// P += (1-alpha_t-1) * pdf_normal(b, mean[state], sd[state]);
	// b = BAF value

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
			// Get the normalized PDF value
			double normalized_pdf = pdf_normalization(b, mean50_state1, sd50_state1);

			// Update the emission probability
			p += (1 - uf) * normalized_pdf;
			// p += (1 - uf) * pdf_normal(b, mean50_state1, sd50_state1);
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
			// Get the normalized PDF value
			double normalized_pdf_left = pdf_normalization(b, mean0, sd0);
			double normalized_pdf_right = pdf_normalization(b, 1 - mean0, sd0);

			// double p_update = p;
			// p_update += ((1 - uf) * (1 - pfb) * normalized_pdf_left) + ((1 - uf) * pfb * normalized_pdf_right);

			// // If p is greater than 1, throw an error and show the equation
			// // values
			// if (p_update > 1)
			// {
			// 	std::cerr << "ERROR [1] = " << p_update << std::endl;
			// 	std::cerr << "UF = " << uf << std::endl;
			// 	std::cerr << "PFB = " << pfb << std::endl;
			// 	std::cerr << "B = " << b << std::endl;
			// 	std::cerr << "State = " << state << std::endl;
			// 	std::cerr << "Equation = " << (1 - uf) << " * (1 - " << pfb << ") * " << normalized_pdf_left << " + " << (1 - uf) << " * " << pfb << " * " << normalized_pdf_right << std::endl;
			// 	std::cerr << "First part = " << (1 - uf) * (1 - pfb) * normalized_pdf_left << std::endl;
			// 	std::cerr << "Second part = " << (1 - uf) * pfb * normalized_pdf_right << std::endl;
			// }

			// Update the emission probability
			p += (1 - uf) * (1 - pfb) * normalized_pdf_left;
			p += (1 - uf) * pfb * normalized_pdf_right;
			// p += (1 - uf) * (1 - pfb) * pdf_normal(b, mean0, sd0);
			// p += (1 - uf) * pfb * pdf_normal(b, 1 - mean0, sd0);
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
			// Get the normalized PDF values
			double normalized_pdf_left = pdf_normalization(b, mean0, sd0);
			double normalized_pdf_middle = pdf_normalization(b, mean50, sd50);
			double normalized_pdf_right = pdf_normalization(b, 1 - mean0, sd0);

			// Update the emission probability
			p += (1 - uf) * (1 - pfb) * (1 - pfb) * normalized_pdf_left;
			p += (1 - uf) * 2 * pfb * (1 - pfb) * normalized_pdf_middle;
			p += (1 - uf) * pfb * pfb * normalized_pdf_right;

			// p += (1 - uf) * (1 - pfb) * (1 - pfb) * pdf_normal(b, mean0, sd0);
			// p += (1 - uf) * 2 * pfb * (1 - pfb) * pdf_normal(b, mean50, sd50);
			// p += (1 - uf) * pfb * pfb * pdf_normal(b, 1 - mean0, sd0);
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
			// Get the normalized PDF values
			double normalized_pdf_left = pdf_normalization(b, mean0, sd0);
			double normalized_pdf_right = pdf_normalization(b, 1 - mean0, sd0);

			// Update the emission probability
			p += (1 - uf) * (1 - pfb) * normalized_pdf_left;
			p += (1 - uf) * pfb * normalized_pdf_right;
			// p += (1 - uf) * (1 - pfb) * pdf_normal(b, mean0, sd0);
			// p += (1 - uf) * pfb * pdf_normal(b, 1 - mean0, sd0);
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
			// Get the normalized PDF values
			double normalized_pdf_1 = pdf_normalization(b, mean0, sd0);
			p += (1 - uf) * (1 - pfb) * (1 - pfb) * (1 - pfb) * normalized_pdf_1;

			double normalized_pdf_2 = pdf_normalization(b, mean33, sd33);
			p += (1 - uf) * 3 * (1 - pfb) * (1 - pfb) * pfb * normalized_pdf_2;

			double normalized_pdf_3 = pdf_normalization(b, 1 - mean33, sd33);
			p += (1 - uf) * 3 * (1 - pfb) * pfb * pfb * normalized_pdf_3;

			double normalized_pdf_4 = pdf_normalization(b, 1 - mean0, sd0);
			p += (1 - uf) * pfb * pfb * pfb * normalized_pdf_4;

			// p += (1 - uf) * (1 - pfb) * (1 - pfb) * (1 - pfb) * pdf_normal(b, mean0, sd0);
			// p += (1 - uf) * 3 * (1 - pfb) * (1 - pfb) * pfb * pdf_normal(b, mean33, sd33);
			// p += (1 - uf) * 3 * (1 - pfb) * pfb * pfb * pdf_normal(b, 1 - mean33, sd33);
			// p += (1 - uf) * pfb * pfb * pfb * pdf_normal(b, 1 - mean0, sd0);
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
			// Get the normalized PDF values
			double normalized_pdf_1 = pdf_normalization(b, mean0, sd0);
			p += (1 - uf) * (1 - pfb) * (1 - pfb) * (1 - pfb) * (1 - pfb) * normalized_pdf_1;

			double normalized_pdf_2 = pdf_normalization(b, mean25, sd25);
			p += (1 - uf) * 4 * (1 - pfb) * (1 - pfb) * (1 - pfb) * pfb * normalized_pdf_2;

			double normalized_pdf_3 = pdf_normalization(b, mean50, sd50);
			p += (1 - uf) * 6 * (1 - pfb) * (1 - pfb) * pfb * pfb * normalized_pdf_3;

			double normalized_pdf_4 = pdf_normalization(b, 1 - mean25, sd25);
			p += (1 - uf) * 4 * (1 - pfb) * pfb * pfb * pfb * normalized_pdf_4;

			double normalized_pdf_5 = pdf_normalization(b, 1 - mean0, sd0);
			p += (1 - uf) * pfb * pfb * pfb * pfb * normalized_pdf_5;

			// p += (1 - uf) * (1 - pfb) * (1 - pfb) * (1 - pfb) * (1 - pfb) * pdf_normal(b, mean0, sd0);
			// p += (1 - uf) * 4 * (1 - pfb) * (1 - pfb) * (1 - pfb) * pfb * pdf_normal(b, mean25, sd25);
			// p += (1 - uf) * 6 * (1 - pfb) * (1 - pfb) * pfb * pfb * pdf_normal(b, mean50, sd50);
			// p += (1 - uf) * 4 * (1 - pfb) * pfb * pfb * pfb * pdf_normal(b, 1 - mean25, sd25);
			// p += (1 - uf) * pfb * pfb * pfb * pfb * pdf_normal(b, 1 - mean0, sd0);
		}
	}
	// if (p == 0)  // Prevent divide by zero error
	// {
	// 	p = FLOAT_MINIMUM;
	// }

	// Ensure that p is between FLOAT_MINIMUM and PROB_MAX
	p = std::max(FLOAT_MINIMUM, std::min(PROB_MAX, p));

	// // Check if p is greater than 1
	// if (p > 1)
	// {
	// 	std::cerr << "ERROR [4] = " << p << std::endl;
	// 	std::cerr << "UF = " << uf << std::endl;
	// 	std::cerr << "PFB = " << pfb << std::endl;
	// 	std::cerr << "B = " << b << std::endl;
	// 	std::cerr << "State = " << state << std::endl;
	// }

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

// SV calling with the HMM via the Viterbi algorithm
std::pair<std::vector<int>, double> ViterbiLogNP_CHMM(CHMM hmm, int T, std::vector<double>& O1, std::vector<double>& O2, std::vector<double>& pfb, double **delta, int **psi, std::vector<double>& pprob)
{
	int i, j; /* state indices */
	int t;	  /* time index */

	int snp_count = 0;

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
			A1[i][j] = hmm.A[i][j];
		}
	}

	/* 0. Preprocessing */

	// Pi is the initial probability distribution over states (The probability that the Markov chain will start in state i).
	// Threshold any zero values to avoid calculation issues.
	for (i = 1; i <= hmm.N; i++)
	{
		if (hmm.pi[i] == 0)
			hmm.pi[i] = 1e-9; /*eliminate problems with zero probability*/
		hmm.pi[i] = log(hmm.pi[i]);  // Convert to log probability due to underflow
	}

	// Bi(Ot) is an NxT matrix of emission probabilities (observation likelihoods)
	// expressing the probability of an observation Ot being generated from a
	// state i. Ot is the observation symbol at time t (in this case, the LRR
	// and BAF values).
	// Initialize the emission probability matrix
	biot = dmatrix(1, hmm.N, 1, T);  // Allocate a NxT double matrix (N=6 states)

	// Determine the minimum and maximum values of the LRR and BAF PDFs for
	// normalization.
	// double min_pdf = 0.0;
	// double max_pdf_b1 = pdf_normal(hmm.B1_mean[1], hmm.B1_mean[1], hmm.B1_sd[1]);
	// double max_pdf_b2 = pdf_normal(hmm.B2_mean[1], hmm.B2_mean[1], hmm.B2_sd[1]);

	//std::cout << "[HMM] Running Viterbi algorithm with " << hmm.N << " states and " << T << " probes\n";

	// Loop through each state N
	// Start at 1 because states are 1-based (1-6)
	for (i = 1; i <= hmm.N; i++)
	{

		// Loop through each probe T (0-based)
		for (t = 1; t <= T; t++)
		// for (t = 0; t < T; t++)
		{
			// A regular SNP marker; we use both LRR and BAF information to
			// calculate the joint probability of the marker being in state i

			// Get the state observation likelihood b_j(O_t) of the observation
			// symbol O_t given the current state j
			// B1_uf is the previous alpha (transition probability)
			double O1_val = O1[t];
			double O1_logprob = b1iot(i, hmm.B1_mean, hmm.B1_sd, hmm.B1_uf, O1_val);

			double O2_val = O2[t];
			double O2_logprob = b2iot(i, hmm.B2_mean, hmm.B2_sd, hmm.B2_uf, pfb[t], O2_val);
			// double O2_logprob = b1iot(i, hmm.B2_mean, hmm.B2_sd, hmm.B2_uf, O2_val);

			// Update the emission probability matrix with the joint probability
			// of the marker being in state i at time t (log probability) based
			// on both LRR and BAF values
			biot[i][t] = O1_logprob + O2_logprob;

			// If this value is greater than 0, print an error
			// if (biot[i][t] > 0)
			// {
			// 	std::cerr << "ERROR [1] = " << biot[i][t] << std::endl;
			// 	std::cerr << "O1_logprob = " << O1_logprob << std::endl;
			// 	std::cerr << "O2_logprob = " << O2_logprob << std::endl;
			// }

			// Update the SNP count
			snp_count++;
		}
	}

	// Loop through biot and check for any positive values
	// for (i = 1; i <= hmm.N; i++)
	// {
	// 	for (t = 1; t <= T; t++)
	// 	{
	// 		if (biot[i][t] > 0)
	// 		{
	// 			std::cerr << "ERROR [2] = " << biot[i][t] << std::endl;
	// 			std::cerr << "i = " << i << std::endl;
	// 			std::cerr << "t = " << t << std::endl;
	// 			std::cerr << "T = " << T << std::endl;
	// 		}
	// 	}
	// }

	/*fprintf(stderr, "NOTICE: Encounterd %i snp probes and %i cn probes\n", snp_count, cn_count);*/

	/* 1. Initialization  */

	for (i = 1; i <= hmm.N; i++)
	{
		delta[1][i] = hmm.pi[i] + biot[i][1];  // Initialize the delta matrix (log probability)
		psi[1][i] = 0;  // Initialize the psi matrix (state sequence) to 0 (no state)
		pprob[i] = -VITHUGE;  // Initialize posterior log probabilities for each state to -inf
	}

	/* 2. Recursion */
	for (t = 2; t <= T; t++)
	{

		// Based on the SNP distance, update the transition matrix
		// Not used in this implementation (snpdist is always 1)
		// if (hmm.dist != 1)
		// 	convertHMMTransition(hmm, A1, snpdist[t - 1]); /*t-1 is used because the current val is calculated from previous values*/

		for (j = 1; j <= hmm.N; j++)
		{
			maxval = -VITHUGE;
			maxvalind = 1;
			for (i = 1; i <= hmm.N; i++)
			{
				val = delta[t - 1][i] + log(A1[i][j]);  // Update the delta matrix (log probability)
				if (val > maxval)  // Update the max value
				{
					maxval = val;
					maxvalind = i;
				}
			}

			delta[t][j] = maxval + biot[j][t];  // Update the delta matrix (log probability)

			// If this value is greater than 0, print an error
			// if (delta[t][j] > 0)
			// {
			// 	std::cerr << "ERROR [0] = " << delta[t][j] << std::endl;
			// 	std::cerr << "maxval = " << maxval << std::endl;
			// 	std::cerr << "biot = " << biot[j][t] << std::endl;
			// }

			psi[t][j] = maxvalind;  // Update the psi matrix (state sequence)
		}
	}

	/* 3. Termination */

	// Variable to store the most likely state at time T (maximum a posteriori
	// probability estimate)
	double max_a_posteriori_prob = -VITHUGE;

	// *pprob = -VITHUGE;
	q[T] = 1;
	for (i = 1; i <= hmm.N; i++)
	{
		// Get the log probability of state i at time T
		double current_prob = delta[T][i];

		// Get the posterior log probability of the state
		double prev_prob = pprob[i];

		// If it's greater than the initial probability, update the probability
		if (current_prob > prev_prob)
		{
			// *pprob = prob;
			q[T] = i;  // Set the state at time T to i
			pprob[i] = current_prob;  // Update the posterior log probability for state i

			// Update the maximum a posteriori probability estimate
			max_a_posteriori_prob = current_prob;
		}
	}

	/* 4. Path (state sequence) backtracking */
	for (t = T - 1; t >= 1; t--)
		q[t] = psi[t + 1][q[t + 1]];

	for (i = 1; i <= hmm.N; i++)
	{ /*recover the HMM model as original*/
		hmm.pi[i] = exp(hmm.pi[i]);
	}
	free_dmatrix(biot, 1, hmm.N, 1, T);
	free_dmatrix(A1, 1, hmm.N, 1, hmm.N);

	// Return the state sequence and its likelihood
	return std::make_pair(q, max_a_posteriori_prob);
}

CHMM ReadCHMM(const char *filename)
{
	FILE *fp;
	CHMM hmm;
	int i, j, k;

	fp = fopen(filename, "r");
	if (!fp)
		fprintf(stderr, "Error: cannot read from HMM file %s\n", filename);

	if (fscanf(fp, "M=%d\n", &(hmm.M)) == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read M annotation from HMM file");
	if (fscanf(fp, "N=%d\n", &(hmm.N)) == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read N annotation from HMM file");

	if (fscanf(fp, "A:\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read A annotation from HMM file");
	hmm.A = (double **)dmatrix(1, hmm.N, 1, hmm.N);
	for (i = 1; i <= hmm.N; i++)
	{
		for (j = 1; j <= hmm.N; j++)
		{
			if (fscanf(fp, "%lf", &(hmm.A[i][j])) == EOF)
				fprintf(stderr, "khmm::ReadCHMM: cannot read A matrix from HMM file");
		}
		if (fscanf(fp, "\n") == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");
	}

	if (fscanf(fp, "B:\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read B annotation from HMM file");
	hmm.B = (double **)dmatrix(1, hmm.N, 1, hmm.M);
	for (j = 1; j <= hmm.N; j++)
	{
		for (k = 1; k <= hmm.M; k++)
		{
			if (fscanf(fp, "%lf", &(hmm.B[j][k])) == EOF)
				fprintf(stderr, "khmm::ReadCHMM: cannot read B matrix from HMM file");
		}
		if (fscanf(fp, "\n") == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");
	}

	if (fscanf(fp, "pi:\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read PI annotation from HMM file");
	hmm.pi = (double *)dvector(1, hmm.N);
	for (i = 1; i <= hmm.N; i++)
	{
		if (fscanf(fp, "%lf", &(hmm.pi[i])) == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read PI vector from HMM file");
		if (hmm.pi[i] < 1e-6)
			hmm.pi[i] = 1e-6;
	}
	if (fscanf(fp, "\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");

	if (fscanf(fp, "B1_mean:\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read B1_mean annotation from HMM file");
	hmm.B1_mean = (double *)dvector(1, hmm.N);
	for (i = 1; i <= hmm.N; i++)
		if (fscanf(fp, "%lf", &(hmm.B1_mean[i])) == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read B1_mean vector from HMM file");
	if (fscanf(fp, "\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");

	if (fscanf(fp, "B1_sd:\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read B1_sd annotation from HMM file");
	hmm.B1_sd = (double *)dvector(1, hmm.N);
	for (i = 1; i <= hmm.N; i++)
		if (fscanf(fp, "%lf", &(hmm.B1_sd[i])) == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read B1_sd from HMM file");
	if (fscanf(fp, "\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");

	if (fscanf(fp, "B1_uf:\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read B1_uf annotation from HMM file");
	if (fscanf(fp, "%lf", &(hmm.B1_uf)) == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read B1_uf from HMM file");
	if (fscanf(fp, "\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");

	if (fscanf(fp, "B2_mean:\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read B2_mean annotation from HMM file");
	hmm.B2_mean = (double *)dvector(1, 5);
	for (i = 1; i <= 5; i++)
		if (fscanf(fp, "%lf", &(hmm.B2_mean[i])) == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read B2_mean from HMM file");
	if (fscanf(fp, "\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");

	if (fscanf(fp, "B2_sd:\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read B2_sd annotation from HMM file");
	hmm.B2_sd = (double *)dvector(1, 5);
	for (i = 1; i <= 5; i++)
		if (fscanf(fp, "%lf", &(hmm.B2_sd[i])) == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read B2_sd from HMM file");
	if (fscanf(fp, "\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");

	if (fscanf(fp, "B2_uf:\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read B2_uf annotation from HMM file");
	if (fscanf(fp, "%lf", &(hmm.B2_uf)) == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read B2_uf from HMM file");
	if (fscanf(fp, "\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");

	if (fscanf(fp, "B3_mean:\n") != EOF)
	{
		hmm.NP_flag = 1;
		hmm.B3_mean = (double *)dvector(1, hmm.N);
		for (i = 1; i <= hmm.N; i++)
			if (fscanf(fp, "%lf", &(hmm.B3_mean[i])) == EOF)
				fprintf(stderr, "khmm::ReadCHMM: cannot read B3_mean from HMM file");
		if (fscanf(fp, "\n") == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");
		if (fscanf(fp, "B3_sd:\n") == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read B3_sd annotation from HMM file");
		hmm.B3_sd = (double *)dvector(1, hmm.N);
		for (i = 1; i <= hmm.N; i++)
			if (fscanf(fp, "%lf", &(hmm.B3_sd[i])) == EOF)
				fprintf(stderr, "khmm::ReadCHMM: cannot read B3_sd from HMM file");
		if (fscanf(fp, "\n") == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");
		if (fscanf(fp, "B3_uf:\n") == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read B3_uf annotation from HMM file");
		if (fscanf(fp, "%lf", &(hmm.B3_uf)) == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read B3_uf from HMM file");
		if (fscanf(fp, "\n") == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");
	}
	else
	{
		hmm.NP_flag = 0;
	}

	if (fscanf(fp, "DIST:\n") != EOF)
	{
		if (fscanf(fp, "%d", &(hmm.dist)) == EOF)
			fprintf(stderr, "khmm:ReadCHMM: cannot read DIST from HMM file");
	}
	else
	{
		// hmm.dist = STATE_CHANGE;
		//  snp_dist is the default distance between two SNPs in the same state
		//  (not used in this implementation)
		//  Set it to 1 to disable the distance model
		hmm.dist = 1;
	}

	fclose(fp);
	return hmm;
}

// void FreeCHMM(CHMM *hmm)
// {
// 	free_dmatrix(hmm.A, 1, hmm.N, 1, hmm.N);
// 	free_dmatrix(hmm.B, 1, hmm.N, 1, hmm.M);
// 	free_dvector(hmm.pi, 1, hmm.N);
// 	free_dvector(hmm.B1_mean, 1, hmm.N);
// 	free_dvector(hmm.B1_sd, 1, hmm.N);
// 	free_dvector(hmm.B2_mean, 1, 5);
// 	free_dvector(hmm.B2_sd, 1, 5);

// 	if (hmm.NP_flag)
// 	{
// 		free_dvector(hmm.B3_mean, 1, hmm.N);
// 		free_dvector(hmm.B3_sd, 1, hmm.N);
// 	}
// }
