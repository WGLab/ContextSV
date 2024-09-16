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

typedef struct {
	int N;			/* number of states;  Q={1,2,...,N} */
	int M; 			/* number of observation symbols; V={1,2,...,M}*/
	double **A;		/* A[1..N][1..N]. a[i][j] is the transition prob
			   	of going from state i at time t to state j
			   	at time t+1 */
	double **B;		/* B[1..N][1..M]. b[j][k] is the probability of
			   	of observing symbol k in state j */
	double *pi;		/* pi[1..N] pi[i] is the initial state distribution. */
	double *B1_mean;	/* B1_mean[1..N] mean of a continuous Gaussian distribution for state 1 through N*/
	double *B1_sd;		/*B1_sd standard deviation of B1 values, which is the same for all states*/
	double B1_uf;		/*B1_uniform_fraction: the contribution of uniform distribution to the finite mixture model */
	double *B2_mean;	/* B2_mean[1..4] is the average of B_allele_freq*/
	double *B2_sd;		/* B2_sd[1..4] is the standard deviation of four B_allele_freq, B2_sd[5] is specially for state1, where B is modelled as a wide normal distribution */
	double B2_uf;		/* B2_uniform_fraction: the fraction of uniform distribution in the finite mixture model */
	
	int NP_flag;		/*flag of 1 and 0 to indicate whether Non-Polymorhpic marker information is contained with HMM file*/
	double *B3_mean;	/* B3_mean[1..N] mean of non-polymorphic probe for state 1 through N*/
	double *B3_sd;		/* B3_sd[1..4] is the standard deviation of B3 values*/
	double B3_uf;		/* B3_uniform_fraction: */
	int dist;		/* new parameter to facilitate CNV calling from resequencing data (2009 April) */
} CHMM;


/************************************
*subroutines for processing continuous HMM models
************************************/

/// Read an HMM from a file
CHMM ReadCHMM (const char *filename);

// /// Free the memory allocated for an HMM
// void FreeCHMM(CHMM *phmm);

/// Run the main HMM algorithm
std::pair<std::vector<int>, double> testVit_CHMM(CHMM hmm, int T, std::vector<double>& O1, std::vector<double>& O2, std::vector<double>& pfb);

/// Viterbi algorithm
std::pair<std::vector<int>, double> ViterbiLogNP_CHMM(CHMM phmm, int T, std::vector<double>& O1, std::vector<double>& O2, std::vector<double>& pfb, double **delta, int **psi, double *pprob);

/// O1 emission probability
double b1iot (int state, double *mean, double *sd, double uf, double o);

/// O2 emission probability
double b2iot (int state, double *mean, double *sd, double uf, double pfb, double b);

/// Return the probability of observing a value in a normal distribution, normalized to a range of [min_pdf, max_pdf]
double pdf_normalization(double obs, double mean, double sd);

#endif // _HMM_H_
