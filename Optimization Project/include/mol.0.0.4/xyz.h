#ifndef _XYZ_H_
#define _XYZ_H_

/**
	\file xyz.h
	This file contains structures and functions
	that are used to initialize the atomgrp
	structure from the xyz file.
*/

/**
	Returns a pointer to an xyz struct,
	read from xyz file at path.
	Exits if file does not exist.
*/
struct atomgrp* read_xyz (const char* path, struct prms* prms);
/**
  reads xyz file with extra sa column.
*/
struct atomgrp* read_oldxyz (const char* path, struct prms* prms);

/**
	Prints an xyz file from the xyz coordinates.
*/
void print_xyz (struct atomgrp* xyz, struct prms* prms);
void fprint_xyz (struct atomgrp* xyz, struct prms* prms, const char* path);
void fprint_oldxyz (struct atomgrp* xyz, struct prms* prms, const char* path);

#endif
