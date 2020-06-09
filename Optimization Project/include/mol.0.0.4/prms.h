#ifndef _PRMS_H_
#define _PRMS_H_

/**
	\file prms.h
	This file contains structures and functions
	that are used to initialize prms from
	the atom parameter file.
*/

struct prms
{
	int natom_types; /**< number of atom types */
	int nsubatom_types; /**< number of subatom types */

	char** atypemaj; /**< major names of atom types */
	char** atypemin; /**< minor names of atom types */
	int* subatoms; // subatoms[atom_typen] returns corresponding subatom type number
	float* rs; /**< array of vdw radii, 0 based, rs[(atom type number)] returns radius of that atom type number */
	float* chrgs; /**< array of charges, 0 based, chs[(atom type number)] returns charge of that atom type number */

	struct atype* atypeh; /**< hash (array) of key (atypemaj atypemin) value (atom type number) pairs */

	int hydrogen_typen; // hydrogen type number

	float r1; /**< r1 for potential. */
	float r2; /**< r2 for potential. */

	int k; /**< number of eigenvalues */
	float* lambdas; /**< pointer to array of eigenvalues */
	float* Xs; /**< pointer to matrix of eigenvectors */
};

struct atype // allow lookup of atypen from atypemaj and atypemin
{
	char* key; /**< "atypemaj atypemin" */
	int val; /**< atom type number */
};
/**
	Current line type being read.
*/
enum ereadstate
{
	VERSION,
	ATOM,
	HYDROGEN,
	RADIUS,
	POTENTIAL,
};

enum ereaderr
{
	ERR_VERSION,
	ERR_ATOM,
	ERR_NAMELEN,
	ERR_HYDROGEN,
	ERR_RADIUS,
	ERR_POTENTIAL,
};

/**
	Defines atoms types, pairwise energies,
	and charges.
*/
struct prms* read_prms (const char* path, const char* bin_version);

/**
  Return the atom type value of the atypemaj, atypemin key
*/
int atype_val (struct prms* prms, const char* atypemaj, const char* atypemin);

/**
	Frees the contents of prms and prms itself.
*/
void free_prms (struct prms* prms);

/**
	Prints an error message based on readerr.
*/
void print_readerr (enum ereaderr readerr, const char* path, char* line);

/**
	Prints the prms struct to stdout
*/
void print_prms (struct prms* prms);

#endif
