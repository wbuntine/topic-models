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
#ifndef _EHASH_H
#define _EHASH_H

#include <stdint.h>

/*
 *   stored as +1 to make 0 be "empty"
 */
typedef struct ehash_s {
	uint32_t *tab;
	uint32_t size;
} ehash_t;

int ehash_init(ehash_t *hp, int size);
void ehash_free(ehash_t *hp);
/*
 *    w = word to find
 *    wind[] = index of words
 */
int32_t ehash_findw(ehash_t *hp, uint32_t w, uint32_t *wind);
void ehash_addw(ehash_t *hp, uint32_t w, uint32_t ind);

#endif
