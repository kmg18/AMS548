#ifndef _IO_H_
#define _IO_H_

/** \file io.h
	This file contains structures and functions
	that are used for reading from and
	writing to files
*/

/**
  determines the efile ext of path
  \return an integer representing the file type
*/
int file_ext (const char* path);

/**
  prints to stderr a list of known file types
*/
void fprintf_stderr_atomgrp_file_exts ();

/**
  reads an atomgrp file based on its ext
*/
struct atomgrp* read_file_atomgrp (const char* path, struct prms* prms);

/**
  prints an atomgrp file based on its ext
*/
void fprint_file_atomgrp (const char* path, struct atomgrp* ag, struct prms* prms);

#endif
