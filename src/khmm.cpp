#include <iostream>

#include "khmm.h"
#include "kc.h"

#define STATE_CHANGE 100000.0 /*this is the expected changes (D value) in the transition matrix*/
#define VITHUGE 100000000000.0
#define FLOAT_MINIMUM 1.175494351e-38; /*this is indeed machine dependent*/
#define DELTA 1

/*	This file was re-written from several subroutines from the UMDHMM package by Tapas Kanungo (Date: 15 December 1997), which has excellent framework of the implementation of Forward-Backward, Viterbi, and Baum-Welch algorithms.
	The original UMDHMM package was downloaded from http://www.kanungo.com/software/software.html. The citation for the UMDHMM program is "UMDHMM: Hidden Markov Model Toolkit," in "Extended Finite State Models of Language," A. Kornai (editor), Cambridge University Press, 1999."
	The basic framework (including variable name, subroutine name) is highly similar to the original UMDHMM package, but the actual implementation is completely different as no "discrete symbol emission" is used in PennCNV.
*/

// Entry point
std::vector<int> testVit_CHMM(CHMM hmm, int T, double *O1, double *O2, double *pfb, int *snpdist, double *plogproba)
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

	// int *q;			/* Array for the state sequence q[1..T] */
	double **delta; // Matrix
	int **psi;		// Matrix

	// These 3 variables will be updated with the HMM results
	// q = ivector(1, T);                // Allocate an int vector for each
	// probe
	delta = dmatrix(1, T, 1, hmm.N); // Allocate a TxN  double matrix (N=6 states)
	psi = imatrix(1, T, 1, hmm.N);	 // Allocate a TxN  int matrix (N=6 states)

	// Allocate pfb and plogproba
	pfb = dvector(1, T);
	plogproba = dvector(1, hmm.N);

	// Initialize pfb and plogproba
	for (int i = 1; i <= T; i++)
	{
		pfb[i] = 1;
	}

	for (int i = 1; i <= hmm.N; i++)
	{
		plogproba[i] = -VITHUGE;
	}

	// Run the HMM
	std::vector<int> q;  // State sequence
	q = ViterbiLogNP_CHMM(&hmm, T, O1, O2, pfb, snpdist, delta, psi, plogproba);

	// Free the variables
	// free_ivector(q, 1, T);
	free_imatrix(psi, 1, T, 1, hmm.N);
	free_dmatrix(delta, 1, T, 1, hmm.N);

	// Pop the first element of q, which is always 0 (Done this way for 1-based indexing)
	q.erase(q.begin());

	// Return the state sequence
	return q;
}

// Emission probability: Calculate the state observationn likelihood b_j(O_t) of the observation symbol O_t given the current state j
double b1iot(int state, double *mean, double *sd, double uf, double o)
{
	// UF = previous alpha
	double p = 0;
	p = uf;

	// PDF normal is the transition probability distrubution a_ij (initialized as pi_n) from state i to j
	// P += (1-alpha_t-1) *
	p += (1 - uf) * pdf_normal(o, mean[state], sd[state]);

	// Prevent divide by zero error
	if (p == 0)
		p = FLOAT_MINIMUM;

	// Return the log probability
	return log(p);
}

double b2iot(int state, double *mean, double *sd, double uf, double pfb, double b)
{
	double p = 0;
	double mean0 = mean[1];
	double mean25 = mean[2];
	double mean33 = mean[3];
	double mean50 = mean[4];
	double mean50_state1 = mean[5];
	double sd0 = sd[1];
	double sd25 = sd[2];
	double sd33 = sd[3];
	double sd50 = sd[4];
	double sd50_state1 = sd[5];

	p = uf;
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
			p += (1 - uf) * pdf_normal(b, mean50_state1, sd50_state1);
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
			p += (1 - uf) * (1 - pfb) * pdf_normal(b, mean0, sd0);
			p += (1 - uf) * pfb * pdf_normal(b, 1 - mean0, sd0);
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
			p += (1 - uf) * (1 - pfb) * (1 - pfb) * pdf_normal(b, mean0, sd0);
			p += (1 - uf) * 2 * pfb * (1 - pfb) * pdf_normal(b, mean50, sd50);
			p += (1 - uf) * pfb * pfb * pdf_normal(b, 1 - mean0, sd0);
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
			p += (1 - uf) * (1 - pfb) * pdf_normal(b, mean0, sd0);
			p += (1 - uf) * pfb * pdf_normal(b, 1 - mean0, sd0);
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
			p += (1 - uf) * (1 - pfb) * (1 - pfb) * (1 - pfb) * pdf_normal(b, mean0, sd0);
			p += (1 - uf) * 3 * (1 - pfb) * (1 - pfb) * pfb * pdf_normal(b, mean33, sd33);
			p += (1 - uf) * 3 * (1 - pfb) * pfb * pfb * pdf_normal(b, 1 - mean33, sd33);
			p += (1 - uf) * pfb * pfb * pfb * pdf_normal(b, 1 - mean0, sd0);
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
			p += (1 - uf) * (1 - pfb) * (1 - pfb) * (1 - pfb) * (1 - pfb) * pdf_normal(b, mean0, sd0);
			p += (1 - uf) * 4 * (1 - pfb) * (1 - pfb) * (1 - pfb) * pfb * pdf_normal(b, mean25, sd25);
			p += (1 - uf) * 6 * (1 - pfb) * (1 - pfb) * pfb * pfb * pdf_normal(b, mean50, sd50);
			p += (1 - uf) * 4 * (1 - pfb) * pfb * pfb * pfb * pdf_normal(b, 1 - mean25, sd25);
			p += (1 - uf) * pfb * pfb * pfb * pfb * pdf_normal(b, 1 - mean0, sd0);
		}
	}
	if (p == 0)
		p = FLOAT_MINIMUM;
	return log(p);
}

// SV calling with the HMM via the Viterbi algorithm
std::vector<int> ViterbiLogNP_CHMM(CHMM *phmm, int T, double *O1, double *O2, double *pfb, int *snpdist, double **delta, int **psi, double *pprob)
{
	int i, j; /* state indices */
	int t;	  /* time index */

	int snp_count = 0;
	int cn_count = 0;

	int maxvalind;
	double maxval, val;
	double **biot;
	double **A1;

	std::vector<int> q(T + 1);  // state sequence (+1 for 1-based indexing)

	// A1 is the NxN transition probability matrix from state i to j.
	// Initialize with initial values from the provided HMM file
	A1 = dmatrix(1, phmm->N, 1, phmm->N); /*initialize A1 matrix*/
	for (i = 1; i <= phmm->N; i++)
	{
		for (j = 1; j <= phmm->N; j++)
		{
			A1[i][j] = phmm->A[i][j];
		}
	}

	/* 0. Preprocessing */

	// Pi is the initial probability distribution over states (The probability that the Markov chain will start in state i).
	// Threshold any zero values to avoid calculation issues.
	for (i = 1; i <= phmm->N; i++)
	{
		if (phmm->pi[i] == 0)
			phmm->pi[i] = 1e-9; /*eliminate problems with zero probability*/
		phmm->pi[i] = log(phmm->pi[i]);
	}

	// Bi(Ot) is an NxT matrix of emission probabilities (observation likelihoods)
	// expressing the probability of an observation Ot being generated from a state i
	biot = dmatrix(1, phmm->N, 1, T);

	// Loop through each state N
	for (i = 1; i <= phmm->N; i++)
	{

		// Loop through each probe T
		for (t = 1; t <= T; t++)
		{
			// If it's non-polymorphic (BAF > 1)
			if (O2[t] > 1)
			{ /*normally BAF>=0 and BAF<=1; we use BAF>1 to indicate a Non-Polymorphic marker*/
				if (!phmm->NP_flag)
					std::cerr << "FATAL ERROR: CN probe detected but HMM model does not contain parameters for them\n";

				// Update emissions based on LRR
				biot[i][t] = b1iot(i, phmm->B3_mean, phmm->B3_sd, phmm->B3_uf, O1[t]);
				cn_count++;

				// If it's a regular SNP
			}
			else
			{ /*a regular SNP marker; we use both LRR and BAF information to calculate logProb for the marker*/
				// Update emissions based on LRR (O1)
				biot[i][t] = b1iot(i, phmm->B1_mean, phmm->B1_sd, phmm->B1_uf, O1[t]);

				// Update emissions based on BAF (O2)
				// *Set pfb to all 1's to ignore BAF population frequency in
				// this implementation
				biot[i][t] += b2iot(i, phmm->B2_mean, phmm->B2_sd, phmm->B2_uf, pfb[t], O2[t]);
				snp_count++;
			}
		}
	}
	/*fprintf(stderr, "NOTICE: Encounterd %i snp probes and %i cn probes\n", snp_count, cn_count);*/

	/* 1. Initialization  */

	for (i = 1; i <= phmm->N; i++)
	{
		delta[1][i] = phmm->pi[i] + biot[i][1];
		psi[1][i] = 0;
		pprob[1] = -VITHUGE;  // Initialize log probabilities for each state to -inf
	}

	/* 2. Recursion */

	for (t = 2; t <= T; t++)
	{

		// Based on the SNP distance, update the transition matrix
		// Not used in this implementation (snpdist is always 1)
		// if (phmm->dist != 1)
		// 	convertHMMTransition(phmm, A1, snpdist[t - 1]); /*t-1 is used because the current val is calculated from previous values*/

		for (j = 1; j <= phmm->N; j++)
		{
			maxval = -VITHUGE;
			maxvalind = 1;
			for (i = 1; i <= phmm->N; i++)
			{
				val = delta[t - 1][i] + log(A1[i][j]);
				if (val > maxval)
				{
					maxval = val;
					maxvalind = i;
				}
			}

			delta[t][j] = maxval + biot[j][t];
			psi[t][j] = maxvalind;
		}
	}

	/* 3. Termination */

	// *pprob = -VITHUGE;
	q[T] = 1;
	for (i = 1; i <= phmm->N; i++)
	{
		// Get the log probability of state i at time T
		double current_prob = delta[T][i];

		// Get the current log probability of the state
		double prev_prob = pprob[i];

		// If it's greater than the initial probability, update the probability
		if (current_prob > prev_prob)
		{
			// *pprob = prob;
			q[T] = i;  // Set the state at time T to i
			pprob[i] = current_prob;  // Update the log probability of state i
		}
	}

	/* 4. Path (state sequence) backtracking */

	for (t = T - 1; t >= 1; t--)
		q[t] = psi[t + 1][q[t + 1]];

	for (i = 1; i <= phmm->N; i++)
	{ /*recover the HMM model as original*/
		phmm->pi[i] = exp(phmm->pi[i]);
	}
	free_dmatrix(biot, 1, phmm->N, 1, T);
	free_dmatrix(A1, 1, phmm->N, 1, phmm->N);

	// Return the state sequence
	return q;
}

CHMM ReadCHMM(const char *filename)
{
	FILE *fp;
	CHMM hmm;
	CHMM *phmm;
	int i, j, k;

	phmm = &hmm;
	fp = fopen(filename, "r");
	if (!fp)
		fprintf(stderr, "Error: cannot read from HMM file %s\n", filename);

	if (fscanf(fp, "M=%d\n", &(phmm->M)) == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read M annotation from HMM file");
	if (fscanf(fp, "N=%d\n", &(phmm->N)) == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read N annotation from HMM file");

	if (fscanf(fp, "A:\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read A annotation from HMM file");
	phmm->A = (double **)dmatrix(1, phmm->N, 1, phmm->N);
	for (i = 1; i <= phmm->N; i++)
	{
		for (j = 1; j <= phmm->N; j++)
		{
			if (fscanf(fp, "%lf", &(phmm->A[i][j])) == EOF)
				fprintf(stderr, "khmm::ReadCHMM: cannot read A matrix from HMM file");
		}
		if (fscanf(fp, "\n") == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");
	}

	if (fscanf(fp, "B:\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read B annotation from HMM file");
	phmm->B = (double **)dmatrix(1, phmm->N, 1, phmm->M);
	for (j = 1; j <= phmm->N; j++)
	{
		for (k = 1; k <= phmm->M; k++)
		{
			if (fscanf(fp, "%lf", &(phmm->B[j][k])) == EOF)
				fprintf(stderr, "khmm::ReadCHMM: cannot read B matrix from HMM file");
		}
		if (fscanf(fp, "\n") == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");
	}

	if (fscanf(fp, "pi:\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read PI annotation from HMM file");
	phmm->pi = (double *)dvector(1, phmm->N);
	for (i = 1; i <= phmm->N; i++)
	{
		if (fscanf(fp, "%lf", &(phmm->pi[i])) == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read PI vector from HMM file");
		if (phmm->pi[i] < 1e-6)
			phmm->pi[i] = 1e-6;
	}
	if (fscanf(fp, "\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");

	if (fscanf(fp, "B1_mean:\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read B1_mean annotation from HMM file");
	phmm->B1_mean = (double *)dvector(1, phmm->N);
	for (i = 1; i <= phmm->N; i++)
		if (fscanf(fp, "%lf", &(phmm->B1_mean[i])) == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read B1_mean vector from HMM file");
	if (fscanf(fp, "\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");

	if (fscanf(fp, "B1_sd:\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read B1_sd annotation from HMM file");
	phmm->B1_sd = (double *)dvector(1, phmm->N);
	for (i = 1; i <= phmm->N; i++)
		if (fscanf(fp, "%lf", &(phmm->B1_sd[i])) == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read B1_sd from HMM file");
	if (fscanf(fp, "\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");

	if (fscanf(fp, "B1_uf:\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read B1_uf annotation from HMM file");
	if (fscanf(fp, "%lf", &(phmm->B1_uf)) == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read B1_uf from HMM file");
	if (fscanf(fp, "\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");

	if (fscanf(fp, "B2_mean:\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read B2_mean annotation from HMM file");
	phmm->B2_mean = (double *)dvector(1, 5);
	for (i = 1; i <= 5; i++)
		if (fscanf(fp, "%lf", &(phmm->B2_mean[i])) == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read B2_mean from HMM file");
	if (fscanf(fp, "\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");

	if (fscanf(fp, "B2_sd:\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read B2_sd annotation from HMM file");
	phmm->B2_sd = (double *)dvector(1, 5);
	for (i = 1; i <= 5; i++)
		if (fscanf(fp, "%lf", &(phmm->B2_sd[i])) == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read B2_sd from HMM file");
	if (fscanf(fp, "\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");

	if (fscanf(fp, "B2_uf:\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read B2_uf annotation from HMM file");
	if (fscanf(fp, "%lf", &(phmm->B2_uf)) == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read B2_uf from HMM file");
	if (fscanf(fp, "\n") == EOF)
		fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");

	if (fscanf(fp, "B3_mean:\n") != EOF)
	{
		phmm->NP_flag = 1;
		phmm->B3_mean = (double *)dvector(1, phmm->N);
		for (i = 1; i <= phmm->N; i++)
			if (fscanf(fp, "%lf", &(phmm->B3_mean[i])) == EOF)
				fprintf(stderr, "khmm::ReadCHMM: cannot read B3_mean from HMM file");
		if (fscanf(fp, "\n") == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");
		if (fscanf(fp, "B3_sd:\n") == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read B3_sd annotation from HMM file");
		phmm->B3_sd = (double *)dvector(1, phmm->N);
		for (i = 1; i <= phmm->N; i++)
			if (fscanf(fp, "%lf", &(phmm->B3_sd[i])) == EOF)
				fprintf(stderr, "khmm::ReadCHMM: cannot read B3_sd from HMM file");
		if (fscanf(fp, "\n") == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");
		if (fscanf(fp, "B3_uf:\n") == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read B3_uf annotation from HMM file");
		if (fscanf(fp, "%lf", &(phmm->B3_uf)) == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read B3_uf from HMM file");
		if (fscanf(fp, "\n") == EOF)
			fprintf(stderr, "khmm::ReadCHMM: cannot read return character from HMM file");
	}
	else
	{
		phmm->NP_flag = 0;
	}

	if (fscanf(fp, "DIST:\n") != EOF)
	{
		if (fscanf(fp, "%d", &(phmm->dist)) == EOF)
			fprintf(stderr, "khmm:ReadCHMM: cannot read DIST from HMM file");
	}
	else
	{
		// phmm->dist = STATE_CHANGE;
		//  snp_dist is the default distance between two SNPs in the same state
		//  (not used in this implementation)
		//  Set it to 1 to disable the distance model
		phmm->dist = 1;
	}

	fclose(fp);
	return hmm;
}

void FreeCHMM(CHMM *phmm)
{
	free_dmatrix(phmm->A, 1, phmm->N, 1, phmm->N);
	free_dmatrix(phmm->B, 1, phmm->N, 1, phmm->M);
	free_dvector(phmm->pi, 1, phmm->N);
	free_dvector(phmm->B1_mean, 1, phmm->N);
	free_dvector(phmm->B1_sd, 1, phmm->N);
	free_dvector(phmm->B2_mean, 1, 5);
	free_dvector(phmm->B2_sd, 1, 5);

	if (phmm->NP_flag)
	{
		free_dvector(phmm->B3_mean, 1, phmm->N);
		free_dvector(phmm->B3_sd, 1, phmm->N);
	}
}

