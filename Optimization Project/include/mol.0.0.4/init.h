#ifndef _INIT_H_
#define _INIT_H_

/**
	\file init.h
	This file contains functions
	that are used to initialize the
	parameter file locations.
*/

char* main_dir ();

char* prms_dir ();

char* atom_prm_file (const char* atom_prm);
char* coeffs_prm_file (const char* coeffs_prm);
char* rots_prm_file (const char* rots_prm);

#endif
