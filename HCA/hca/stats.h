/*
 * Various data structures for statistics
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

#ifndef __STATS_H
#define __STATS_H

#include "hca.h"
#include "data.h"
#include "misi.h"

/*
 *  all statistics used in the training algorithm
 */
typedef struct D_stats_s {
  /*
   *    this is only used sometimes ('-r phi', '-l phi,...')
   *    and stores the estimated phi matrix
   */
  float **phi;
  float *alpha;
  /*
   *  Basic topic data for simplest LDA model
   *     we assume maximum number of words in single document
   *     fits in 15 bit, use top bit for indicator
   */
  uint16_t *z;
  uint16_t **Ndt;  // number of words in doc for topic
  uint16_t *NdT;   // \sum_t Ndt[d][t]
  uint32_t **Nwt;  // number of words w in topic t
  uint32_t *NWt;    //  \sum_w Nwt[w][t]
  /*
   *  Table counts and statistics thereof for PYP docXtopic models
   */
  uint16_t **Tdt;      // table counts for topics in doc
  uint32_t *TDt, TDT;  // TDt[t] = \sum_d Tdt[d][t] ;  TDT = \sum_t TDt[t]
  uint32_t *Tlife;     //  lifetime in cycles for this
  uint32_t TDTnz;      // number of non-zero TDt[t]
  /*
   *   Table counts and statistics thereof for PYP word models
   *       basic topic by word table counts and its various totals
   *       roughly matches **Tdt
   */
  uint16_t **Twt;   // table counts for words in topic
  uint32_t *TwT, TWT;  // TwT[w] = \sum_t Twt[w][t]; TWT = \sum_w TwT[w]
  uint32_t  TWTnz;  //  = \sum_w 1(TwT[w]>0)
  uint32_t *TWt;  //  TWt[t] = \sum_w Twt[w][t]; 
} D_stats_t;

extern D_stats_t ddS;
extern D_DMi_t ddM;

#define M_multi(l)  misi_multi(&ddM,l)

double gibbs_lda(enum GibbsType fix, int Tmax, int doc, int words, float *p, D_MiSi_t *Dd, int incremental);

/*
 *    steps inside Gibbs to add/remove effects of one word on all stats
 */
void update_topic(int i, int did, int wid, int t, int mi, int *Td_, 
		  D_MiSi_t *Dd, float ttip, float wtip, float dtip);
int remove_topic(int i, int did, int wid, int t, int mi, int *Td_, 
		 D_MiSi_t *Dd, int incremental);

/*
 *    allocation, deallocation, read/write on ddS.z[]
 */
void hca_free();
void hca_alloc();
void hca_rand_z(int Tinit, int firstdoc, int lastdoc);
void hca_read_z(char *resstem, int firstdoc, int lastdoc);
void hca_reset_stats(char *resstem, int restart, int zero,
		     int firstdoc, int lastdoc);
void hca_write_z(char *resstem);
void hca_report(char *resstem, char *stem, int ITER, int procs,
		enum GibbsType fix, int dopmi, int showlike, int nopar);
void hca_correct_twt();

/*
 *     loads data from Chang and Blei's C++ HDP code
 */
void hca_load_hdp(char *resstem);

/*
 *   probabilities for parts of models:
 *      different depending on whether samping or
 *      estimating
 *   generally, *fact() versions for for sampling, *prob() for
 *      estimating
 */
double betabasewordprob(int j);
double wordfact(int j, int t, float *tip);
double wordprob(int j, int t);
void   tableindicatorprob(int j, int t, double *uone, double *uzero);
double remainderindicatorprob(int t);
double alphabasetopicprob(int t);
double topicfact(int d, int t, int tot, uint16_t *zerod, float *tip);
double topicprob(int d, int t, int Ttot);
double topicnorm(int d);
double docfact(D_MiSi_t *dD, int t, int i, int mi, double pK, float *dip);
double docprob(D_MiSi_t *dD, int t, int i, int mi, double pw);

#include "change.h"

#endif
