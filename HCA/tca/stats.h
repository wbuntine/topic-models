/*
 * Various data structures for statistics
 * Copyright (C) 2014 Jinjing Li and Wray Buntine
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Authors: Jinjing Li <jinjingli@gmail.com>
 *          Wray Buntine (wray.buntine@monash.edu)
 *          
 * CHANGES:
 *          don't have [e][d] as an index, just [d]
 *          for s_evt stuff, always put t index last
 */

#ifndef __STATS_H
#define __STATS_H

#include "tca.h"
#include "data.h"
#include "misi.h"

/*
 *  all statistics used in the training algorithm
 */
typedef struct D_stats_s {
  /*
   *  Basic topic data for simplest LDA model
   *     we assume maximum number of words in single document
   *     fits in 15 bit, use top bit for indicator
   */
  uint16_t *z;
  uint16_t **n_dt;     // number of words in doc for topic
  uint16_t *N_dT;      // \sum_t n_dt[d][t]
  uint32_t ***m_evt;   // number of words w in topic t for epoch e
  uint32_t **M_eVt;    //  \sum_w m_evt[e][w][t]  
  /*
   *   the doc PYP (with ad,bd) table totals;
   *   table counts matching ddD.Mi[] and ddD.MI[]
   */
  uint16_t *Si;     // corr. table count (for multis only), indexed by d
  /*
   *  Table counts and statistics thereof for PYP docXtopic models
   */
  uint16_t **c_dt;       // table counts for topics in doc, matches n_dt
  uint16_t *C_dT;        // = \sum_{t} c_dt[d][t]
  uint32_t **C_eDt;      // = \sum_{d\in e} c_dt[d][t]
  uint32_t *C_e;         // = \sum_t C_eDt[e][t]

  /*
   *  Table counts and statistics thereof for tables of PYP docXtopic models
   */
  uint32_t **cp_et;      // table counts for tables of topics in doc
  uint32_t *Cp_e;        // \sum_t cp_et[e][t]

  /*
   *   Table counts and statistics thereof for PYP word models
   *       basic topic by word table counts and its various totals
   *       roughly matches **c_dt
   */
  uint32_t ***s_evt;   // table counts for words in topic
  uint32_t  S_0_nz;    //  = \sum_w 1(S_0vT[w]>0)
  uint32_t  S_0;       //  = \sum_w S_0vT[w]
  uint32_t  *S_0vT;    //  = \sum_t s_evt[0][w][t]
  uint32_t **S_eVt;    //  = \sum_w s_evt[e][w][t];  
} D_stats_t;


extern D_stats_t ddS;
extern D_DMi_t ddM;

#define M_multi(l)  misi_multi(&ddM,l)

double gibbs_lda(enum GibbsType fix, int doc, int words, float *p, D_MiSi_t *dD);
void update_topic(int i, int did, int wid, int t, int mi, 
                  float dtip, D_MiSi_t *dD);
int remove_topic(int i, int did, int wid, int t, int mi, D_MiSi_t *dD);

void get_probs(double *vp);

void tca_report(char *resstem, char *stem, int ITER, int procs, enum GibbsType fix, int dopmi);

void tca_free();
void tca_alloc();
void tca_rand_z(int Tinit, int firstdoc, int lastdoc);
void tca_read_z(char *resstem, int firstdoc, int lastdoc);
void tca_reset_stats(char *resstem, int restart, int warm);
void tca_write_z(char *resstem);

double poch(double b, double a, int N);

void diag_alloc();

/*
 *   probabilities for parts of models
 */
double word_side_fact(int e, int j, int t);
double word_side_prob(int e, int j, int t);
/*
 *   return count to place table back:
 *     0 = no table back
 *     1 = back to previous time
 *     e = back to initial epoch
 *   e+1 = back to root
 *
 *   cuts short if ddP.back in order, so
 *   return is >= ddP.back
 */
int word_side_ind ( int e, int v, int t);

double phi0_prob(int v, int t);

double doc_side_fact(int d, int t);
double doc_side_prob(int d, int t);


int doc_side_ind(int d, int t);

double topicnorm(int d);

double docfact(D_MiSi_t *dD, int t, int i, int mi, double pK, float *dip);
double docprob(D_MiSi_t *dD, int t, int i, int mi, double pw);

/*
 *    during Gibbs estimate proportion
 *    of topics for all documents
 */
void tprob_init();
void tprob_null();
void tprob_free();
void tprob_report(char *resstem, double epsilon);

/*
 *  checks
 */
void check_n_dt(int d);
void tca_checkO();
void tca_checkO_A(int t);
void check_m_evt(int e);
void check_cp_et();

#include "change.h"
#endif
