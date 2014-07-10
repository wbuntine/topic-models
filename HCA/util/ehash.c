/*
 * Dirt simple hashing where size is fixed
 * Copyright (C) 2014 Wray Buntine
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include "yap.h"
#include "ehash.h"

/*********************************
 *    dirt simple hashing routine
 */
#define HASHPRIME 7919U
#define REHASHPRIME 7883U

/*
 *   stored as +1 to make 0 be "empty"
 */

int ehash_init(ehash_t *hp, int size) {
  hp->size = size;
  hp->tab = calloc(hp->size,sizeof(*hp->tab));
  if ( !hp->tab )
    return 1;
  return 0;
}
void ehash_free(ehash_t *hp) {
  free(hp->tab);
  hp->tab = NULL;
  hp->size = 0;
}

/*
 *    w = word to find
 *    wind[] = index of words
 */
int32_t ehash_findw(ehash_t *hp, uint32_t w, uint32_t *wind) {
  uint32_t key = ( w*HASHPRIME ) % hp->size;
  if ( hp->tab[key]==0 ) 
    return UINT32_MAX;
  while ( wind[hp->tab[key]-1] != w ) {
    key = (key+REHASHPRIME) % hp->size;
    if ( hp->tab[key]==0 ) 
      return UINT32_MAX;
  }
  return hp->tab[key]-1;
}
void ehash_addw(ehash_t *hp, uint32_t w, uint32_t ind) {
  uint32_t key = ( w*HASHPRIME ) % hp->size;
  while ( hp->tab[key] != 0 ) {
    key = (key+REHASHPRIME) % hp->size;
  }
  hp->tab[key] = ind+1;
  return;
}

