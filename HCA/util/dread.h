/*
 * Various routines for input data
 * Copyright (C) 2010-2013 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *     
 */

#ifndef __DREAD_H
#define __DREAD_H


enum dataType { WitDit, Docword, LdaC, TxtBag, SeqTxtBag };

typedef struct D_bag_s {
  int W;     // number of unique words
  int D;     // number of docs
  int N;     // number of words in corpus (length of vectors)
  uint32_t   *w, *d;
} D_bag_t;

D_bag_t *data_read(char *stem, enum dataType data);
void data_shrink(D_bag_t *dbp, int size);
void data_append(D_bag_t *dbp, D_bag_t *dbp2);
char *data_name(char *stem, enum dataType data);

#endif
