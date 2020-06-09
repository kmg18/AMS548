#ifndef _PDB_H_
#define _PDB_H_

/**
  read a pdb file
*/
struct atomgrp* read_pdb (const char* path, struct prms* prms);

/**
  print a pdb file
*/
void fprint_pdb (struct atomgrp* ag, struct prms* prms, const char* path);

#endif
