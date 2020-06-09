#ifndef _MYHELPERS_H_
#define _MYHELPERS_H_

/** \file myhelpers.h
	This file contains wrapper functions
	and supporting string functions.
*/

/**
	This is a wrapper to malloc that will print
	an error message if malloc returns NULL.
	\param size allocate this many bytes of memory
	\return a pointer to the allocated memory
*/
void* mymalloc(size_t size);

/**
	This is a wrapper to calloc that will print
	an error message if malloc returns NULL.
	\param nmemb size of array
	\param size size of each element in the array
	\return a pointer to the allocated memory
*/
void* mycalloc(size_t nmemb, size_t size);

/**
	This is a wrapper to realloc that will print
	an error message accordingly.
	\param ptr pointer to the original memory
	\param size allocate this many bytes of memory
	\return a pointer to the allocated memory
*/
void* myrealloc(void* ptr, size_t size);

void myexit (int status);

/**
	This is a wrapper to fopen that will print
	an error message and exit if fopen returns NULL.
*/
FILE* myfopen (const char* path, const char* mode);

/**
	This is a wrapper to fclose that will print
	an error message and exit if an error occurs.
*/
void myfclose (FILE* fp);

/**
	Returns 1 if line consists of whitespace
	characters, else 0;
*/
int iswhiteline (char* line);

#endif
