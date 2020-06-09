#ifndef _COEFFS_H_
#define _COEFFS_H_

/**
	\file coeffs.h
	This file contains structures and functions
	that are used to read in the desired docking
	potential coefficients from the file coeffs.prm.
*/

struct coeff
{
	float rcsw; /**< repulsion correlations weight */
	float nccsw; /**< number of contacts correlations weight */
	float eleccsw; /**< electrostatics correlations weight */
	float borncsw; /**< born electrostatics correlations weight */
	float ecsw; /**< eigen correlations weight */
	float rt; /**< repulsion threshold (max repulsion value allowed in correlation) */
	float nct; /**< ncontacts threshold (min number of contacts allowed in correlation) */
//	float ltbornt; /**< less than born threshold */
//	float gtet; /**< greater than eigen threshold */
};

struct coeffset
{
	int ncoeffs;
	struct coeff* coeffs;
};

struct coeffset* read_coeffs (const char* path);
void free_coeffset (struct coeffset* coeffset);
void print_coeffset (struct coeffset* coeffset);

#endif
