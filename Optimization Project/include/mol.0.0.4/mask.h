#ifndef _MASK_H_
#define _MASK_H_

/**
  mask all atoms in target that are within maskr of any
  atom in mask.
*/
void mask_atomgrp (struct atomgrp* target, struct atomgrp* mask, float maskr);

#endif
