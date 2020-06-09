#ifndef _PROTEIN_H_
#define _PROTEIN_H_

/**
  extracts the residues from ag (as separated by nitrogens)
  and stores the result in an array of atomgrps
*/
struct atomgrp** extract_nitrogen_residues (struct atomgrp* ag, struct prms* prms);

/**
  printf the contents of ress to stdout
*/
void print_residues (struct atomgrp** ress, struct prms* prms);

/**
  free the memory allocated to the array of atomgrps
*/
void free_atomgrps (struct atomgrp** ress);

#endif
