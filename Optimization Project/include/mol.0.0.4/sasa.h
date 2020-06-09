#ifndef _SASA_H_
#define _SASA_H_

void mark_sasa (struct atomgrp* ag, int* sasas);

int* read_sasa (const char* path);

int numbersasas (int* sasas);

#endif
