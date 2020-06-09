#ifndef _MS_H_
#define _MS_H_

/**
  read a marksur pdb file
*/
struct atomgrp* read_ms (const char* path, struct prms* prms);

/**
  print a ms file
*/
void fprint_ms (struct atomgrp* ag, struct prms* prms, const char* path);

#endif
