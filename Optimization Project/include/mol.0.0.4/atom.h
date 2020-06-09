#ifndef _ATOM_H_
#define _ATOM_H_

struct atom
{
	int atom_typen; /**< atom type number */
	int sa; /**< solvent accessible: 1 => solvent accessible, 0 => !1, -1 => undefined */
	float X; /**< X coordinate of atom */
	float Y; /**< Y coordinate of atom */
	float Z; /**< Z coordinate of atom */

	int bonds[4]; /**< list of bonded atoms. if bonds[x] == -1 => bonds[x] doesn't exist */
};

/**
	Holds list of atoms and
	the number of atoms in
	the list.
*/
struct atomgrp
{
	int natoms; /**< number of atoms in the group */
	struct atom* atoms; /**< pointer to the array of atoms in the group */
};

/**
	Frees all atoms in the atomgrp and the atomgrp itself.
*/
void free_atomgrp (struct atomgrp* ag);

/**
	Creates a copy of srcag and returns it.
*/
struct atomgrp* copy_atomgrp (struct atomgrp* srcag);

/**
  Copies contents of atom struct src to dest.
  The pointers must be to previously allocated memory.
*/
void copy_atom (struct atom* src, struct atom* dest);

/**
  extract all atoms of type and return them in atomgrp
*/
struct atomgrp* extract_type (struct atomgrp* ag, const char* type, struct prms* prms);

/**
  remove all atoms of type and return the remaining atoms in atomgrp
*/
struct atomgrp* rm_type (struct atomgrp* ag, const char* type, struct prms* prms);

struct atomgrp* exrm_type (struct atomgrp* ag, const char* type, struct prms* prms, int direction);

/**
  make all ag atoms sa
*/
void full_sa (struct atomgrp* ag);

/**
  join ags into a single atom group ag and return ag
*/
struct atomgrp* join_atomgrps (struct atomgrp** ags);
/**
  join ag1 & ag2 into a single atom group ag and return ag
*/
struct atomgrp* join_2atomgrps (struct atomgrp* ag1, struct atomgrp* ag2);

/**
	Prints the contents of the atomgrp struct to stdout
*/
void print_atomgrp (struct atomgrp* ag, struct prms* prms);

/**
	Prints to stderr the contents of the atom struct
*/
void fprintf_stderr_atom (struct atom* atom, struct prms* prms);

/**
  check the properties of atomgrp for consistency with prms
*/
void check_atomgrp (struct atomgrp* ag, struct prms* prms);

#endif
