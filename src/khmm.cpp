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

	// Free the variables
	free_imatrix(psi, 1, T, 1, hmm.N);
	free_dmatrix(delta, 1, T, 1, hmm.N);

	// Pop the first element of q, which is always 0 (Done this way for 1-based
	// indexing)
	state_sequence.first.erase(state_sequence.first.begin());

	return state_sequence;
}

// Emission probability: Calculate the state observationn likelihood b_j(O_t) of
// the observation symbol O_t given the current state j
// O_t is the LRR value
double b1iot(int state, double *mean, double *sd, double uf, double o)
{
	// Get the normalized PDF value
	double normalized_pdf = pdf_normalization(o, mean[state], sd[state]);

	// Update the emission probability
	double p = uf + ((1 - uf) * normalized_pdf);

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

			// Update the emission probability
			p += (1 - uf) * (1 - pfb) * normalized_pdf_left;
			p += (1 - uf) * pfb * normalized_pdf_right;
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
	biot = dmatrix(1, hmm.N, 1, T);  // Allocate a NxT double matrix (N=6 states)

	// Create a dictionary for each state to store all BAF log probabilities
	std::map<int, std::vector<double>> baf_log_probs;

	// Create a dictionary for each state to store all BAF scaling factors
	std::map<int, std::vector<double>> scaling_factors;

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

// 			// If the BAF value is -1 (no data), set the log probability to 0
// //			double O2_logprob = log(0.5);
// //            double O2_logprob = 0;
// 			if (O2_val == -1) {
// 			    O2_val = 0.5;
// //				O2_logprob = b2iot(i, hmm.B2_mean, hmm.B2_sd, hmm.B2_uf, pfb[t], O2_val);
// 			}


			// If no BAF data (value is -1), only use the LRR data [TEST]
			// double O2_logprob = 0;
			// if (O2_val != -1) {
            // 	O2_logprob = b2iot(i, hmm.B2_mean, hmm.B2_sd, hmm.B2_uf, pfb[t], O2_val);
			// }

			double O2_logprob = b2iot(i, hmm.B2_mean, hmm.B2_sd, hmm.B2_uf, pfb[t], O2_val);

			// Print the scaling factor to determine the contribution of O2 to
			// the emission probability
			// std::cout << "State = " << i << ", O2 = " << O2_val << std::endl;
			// double sum = O1_logprob + O2_logprob;
			// std::cout << "State = " << i <<  ", scaling factor (sum/O1) = "
			// << O1_logprob / sum << std::endl;
			
			// Add the BAF log probability to the dictionary
			baf_log_probs[i].push_back(O2_logprob);

			// Add the scaling factor to the dictionary
			double sum = O1_logprob + O2_logprob;
			scaling_factors[i].push_back(sum / O1_logprob);

			// Update the emission probability matrix with the joint probability
			// of the marker being in state i at time t (log probability) based
			// on both LRR and BAF values
			biot[i][t] = O1_logprob + O2_logprob;

			// Update the SNP count
			snp_count++;
		}
	}

	// Print the mean BAF log probabilities and mean scaling factors for each
	// state
	for (i = 1; i <= hmm.N; i++)
	{
		double sum_baf = 0;
		double sum_scaling = 0;
		for (int j = 0; j < snp_count; j++)
		{
			sum_baf += baf_log_probs[i][j];
			sum_scaling += scaling_factors[i][j];
		}
		double mean_baf = sum_baf / snp_count;
		double mean_scaling = sum_scaling / snp_count;
		std::cout << "State = " << i << ", mean BAF log prob = " << mean_baf
		<< ", mean scaling factor = " << mean_scaling << std::endl;
	}

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
			psi[t][j] = maxvalind;  // Update the psi matrix (state sequence)
		}
	}

	/* 3. Termination */

	// Variable to store the most likely state at time T (maximum a posteriori
	// probability estimate)
	double max_a_posteriori_prob = -VITHUGE;
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
