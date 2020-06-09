#ifndef _ENERGY_H_
#define _ENERGY_H_

/** \file energy.h
	This file contains structures and functions
	that are used for calculating offgrid energies
*/

/**
  find the pairwise potential on the two structures and return it
*/
float pairwise_potential_energy (struct atomgrp* agA, struct atomgrp* agB, struct prms* prms, int only_sab);

/**
  find the coulombic energy the two structures and return it
*/
float coulombic_elec_energy (struct atomgrp* agA, struct atomgrp* agB, struct prms* prms);


/**
  evaluates energy of the 2 interacting proteins. 
*/
float complex_energy (struct atomgrp* agA, struct atomgrp* agB, struct prms* prms);

#endif
