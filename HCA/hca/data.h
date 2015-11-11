/*
 * Various routines for input data
 * Copyright (C) 2010-2014 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@monash.edu)
 *     
 */

#ifndef __DATA_H
#define __DATA_H

#include "hca.h"
#include "dread.h"

/*
 *    Dimensions
 */
typedef struct D_dims_s {
  int W;     // number of unique words
  int C;     // number of classes (got from ddD.c)
  int D, DT; // number of docs
  int DTused; //  DT minus "too small" docs
  int TEST;  // number of docs at end for test
  //  train is 0..D-TEST-1
  int N, NT; // number of words in corpus
  int T; // number of topics
  int DefT;  //   the default topic, might not exist
  char **tokens;  //  optionally loaded
} D_dims_t;

/*
 *    Basic data about documents:
 *    we assume maximum number of words in corpus
 *    fits in 32 bit
 */
typedef struct D_data_s {
  uint32_t *w, *d;
  uint16_t *c;     // class, number in 0,...,C-1
  uint32_t *df;     // document freqency for word
  uint32_t n_df;    // total docs used with df[]
  uint16_t *NdT;    // number of words in doc
  uint16_t NdTmax;  // maximum of NdT[]
  uint32_t *NdTcum; //  cumsum(NdT)
} D_data_t;

extern D_data_t ddD;
extern D_dims_t ddN;

void data_free();
void data_alloc();
int data_docsize();
void data_vocab(char *stem);
void data_class(char *stem);
int data_df(char *stem, uint32_t *df);

void data_report(int ITER, int seed);
void data_checkpoint(char *resstem, char *stem, int ITER);
#endif
