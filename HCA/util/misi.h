/*
 * Module for bursty data preparation and processing
 * Copyright (C) 2013-4 Wray Buntine
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *
 *   dmi_*() routines prepare the collection:
 *           need to record which words occur multiple times in docs
 *   misi_*() routines do per document processing
 *           see dmi_rand() for usage
 */

#ifndef __MISI_H
#define __MISI_H

#include <stdint.h>

/*
 *    used to build global data about word unique usage for a doc
 */
typedef struct D_DMi_s {
  /*
   *  basic dims and data stored/saved from elsewhere
   */
  int T;
  int DT;
  uint32_t *NdTcum;
  uint16_t *z;
  int (*holdout)(int l);   //  call back to say if its holdout
  /*
   *     Only used if doing PDP on docs, i.e., ddP.bdk>0
   *     NB.  keeping stats like this is moderately efficient
   *          and makes sampling ddP.ad and ddP.bdk easier
   */
  uint32_t *multi; //  bit vector saying if w[] has >1 occ. in doc
  uint32_t *multiind; // index to Mi[] and Si[] for non-zero multi[]
  uint16_t *Mi;     // freq. count for multis only
  uint32_t *MI;     // index start of multis for doc (to get singles!)
  int dim_Mi;
  int dim_MiT;       // only include training data
  int dim_multiind;
  int Mi_max;        //  maximum possible Mi[] value, test or train
  int MI_max;        //  maximum number of Mi[] entries for a doc
} D_DMi_t;


/*
 *    used to build local data about topicXword usage for a doc
 */
typedef struct D_MiSi_s {
  /*
   *  where to get globals
   */
  D_DMi_t *pdmi;
  /*
   *    dim 1 of Mik and Sik indicates the word, dim 2 is the topic
   *    the word is indexed like Mi but offset from the smallest
   *    ddD.Mi[] index used by this document ... 
   *    these count and tables (resp.) of this word with topic
   */
  uint16_t **Mik;
  uint16_t **Sik;
  /*
   *  Mi[] = total count for topic in doc,
   *  Si[] = total tables
   */
  uint16_t *Mi;
  uint16_t *Si;
  /*
   *  the smallest ddD.Mi[] index used by this document
   */
  int mi_base;
} D_MiSi_t;

/*
 *   nonzero if this word (index to z[] and d[]) occurs
 *   more than once in its document
 */
#define misi_multi(pdmi,l)  ((pdmi)->multi[(l)/32U] & (1U<<(((unsigned)l)%32U)))

/*
 *  local data handling, 
 *         see misi_rand() code for use
 */
/*  set up a batch/loop over docs  */
void misi_init(D_DMi_t *pdmi, D_MiSi_t *ptr);
void misi_free(D_MiSi_t *ptr);
/*  set up for one document  */
void misi_build(D_MiSi_t *ptr, int d, int notest);
void misi_unbuild(D_MiSi_t *ptr, int d, int notest);
/*  reset stats to zero */
void misi_zero(D_MiSi_t *ptr, int d);
/*  remove/add word from stats */
void misi_decr(D_MiSi_t *ptr, int i, int mi, int t, int w);
int misi_incr(D_MiSi_t *ptr, int i, int mi, int t, int w, float dtip);
/*  word cannot be removed, will violate constraints */
int misi_blocked(D_MiSi_t *ptr, int i, int mi, int t);

/*
 *  global data handling, see elsewhere for use
 */
void dmi_free(D_DMi_t *ptr);
/*
 *      (*ptr) is preallocated memory where the data structure
 *             is built
 *
 *      trainword(l) true if not holdout testing, if if word is not holdout,
 *      use NULL for always true
 *
 *      all other arguments come from the standard data structures
 */
void dmi_init(D_DMi_t *ptr, 
              uint16_t *z, uint32_t *w, uint32_t *NdTcum, //  cumsum(NdT)
              int T, int N, int W, int D, int DT,  // dims
              int (*holdout)(int l));

/*  consistency check */
void dmi_check(D_DMi_t *pdmi, int i);
/*  randomise 'r' assignments for these docs */
void dmi_rand(D_DMi_t *pdmi, int firstdoc, int lastdoc);

#endif
