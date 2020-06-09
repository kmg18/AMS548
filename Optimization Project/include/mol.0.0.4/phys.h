#ifndef _PHYS_H_
#define _PHYS_H_

/** \file phys.h
	This file contains structures and functions
	that are used to
*/


float* moment_of_inertia (struct atomgrp* ag);

void print_moment_of_inertia_matrix (float* moimat);

#endif
