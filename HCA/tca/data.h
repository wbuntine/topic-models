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

#include "tca.h"
#include "dread.h"

/*
 *    Dimensions
 */
typedef struct D_dims_s {
  int W;     // number of unique words
  int E;     // number of epochs
  int D, DT; // number of docs
  int TEST;  // number of docs at end for test
  //  train is 0..D-TEST-1
  int N, NT; // number of words in corpus
  int T;     // number of topics
  int DefT;  //   the default topic, might not exist
  char **tokens;  //  optionally loaded
} D_dims_t;


typedef struct D_data_s {
  /*
   *    Basic data about documents:
   *    we assume maximum number of words in corpus
   *    fits in 32 bit
   */
  uint32_t *w;      //  word,    indexed by n=0,...,N-1
  uint32_t *d;      //  document, indexed by n=0,...,N-1
  uint16_t *e;      //  epoch,   indexed by d=0,...,D-1
  uint16_t *j;      //  doc offset in epoch,   indexed by d=0,...,D-1
  uint32_t *esize;  //  size of epoch,   indexed by e=0,...,E-1
  uint16_t *N_dT;    // number of words in doc
  uint16_t N_dTmax;  // maximum of N_dT[]
  uint32_t *N_dTcum; //  cumsum(N_dT)
} D_data_t;

extern D_data_t ddD;
extern D_dims_t ddN;

/*
 *   return epoch for word in posn i
 */
#define epoch(i)    ddD.e[ddD.d[i]]

void data_free();
void data_alloc();
void data_read_epoch(char *stem);
int data_docsize();
void data_vocab(char *stem);

void data_report(int ITER, int seed);
void data_checkpoint(char *resstem, char *stem, int ITER);
#endif
