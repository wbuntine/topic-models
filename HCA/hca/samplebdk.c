/*
 * Sampling utility for bdk
 * Copyright (C) 2013 Wray Buntine 
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

#include <stdio.h>  
#include <unistd.h>   
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "yap.h"
#include "lgamma.h"
#include "myarms.h"
#include "sample.h"
#include "hca.h"
#include "stats.h"

/*
 *     define to make sample independently for different k
 */
#define SAMPLEBDK_MULTI

/*
 *    docstats[d] points to store for doc d
 *    
 *    stores blocks per nonzero topic
 *       [0] = topic
 *       [1] = total count
 *       [2] = total tables
 *
 *    ending block
 *       [0] = ddN.T+1
 */
static uint16_t **docstats;

/*
 *  globals from hca.h and stats.h
 */
extern double likelihood_bdk();

// #define B_DEBUG

/*
 *     uset<0     ==>   compute for all t
 *     otherwise  ==>   compute for just that t
 */
static double likelihood_bdk_asym(int uset) {
  double *lgbd = dvec(ddN.T);
  double *lgabd = dvec(ddN.T);
  double *lb = dvec(ddN.T);
  int t, i;
  double likelihood = 0;
  if ( uset<0 ) {
    for (t=0; t<ddN.T; t++) {
      lgbd[t] = lgamma(ddP.bdk[t]);
      lgabd[t] = lgamma(ddP.bdk[t]/ddP.ad);
      lb[t] = log(ddP.bdk[t]);
    }
  } else {
    lgbd[uset] = lgamma(ddP.bdk[uset]);
    lgabd[uset] = lgamma(ddP.bdk[uset]/ddP.ad);
    lb[uset] = log(ddP.bdk[uset]);
  }
  for (i=0; i<ddN.DT; i++) {
    uint16_t *store = docstats[i];
    int mem=0;
    if ( store==NULL )
	continue;
    while ( store[mem]<ddN.T ) {
      t = store[mem];
      assert(store[mem+1]>0);
      if ( uset<0 || t==uset ) {
	if ( ddP.ad==0 ) {
	  likelihood += store[mem+2]*lb[t];
	} else {
	  likelihood += 
	    gammadiff(store[mem+2], ddP.bdk[t]/ddP.ad, lgabd[t]);
	}
	likelihood -= gammadiff((int)store[mem+1], ddP.bdk[t], lgbd[t]);
      }
      mem += 3;
    }
    yap_infinite(likelihood);
  }
  if ( uset<0 ) {
    for (t=0; t<ddN.T; t++)
      likelihood += pctl_gammaprior(ddP.bdk[t]);
  } else
    likelihood += pctl_gammaprior(ddP.bdk[uset]);
  free(lgabd);
  free(lgbd);
  free(lb);
  return likelihood;
}

/*
 *    just call the likelihood function
 */
#ifdef SAMPLEBDK_MULTI
static double bdkterms(double b, void *mydata) {
  int t = *((int *)mydata);
  double val;
  ddP.bdk[t] = b;
  val = likelihood_bdk_asym(t);
  myarms_evals++;
#ifdef B_DEBUG
  yap_message("Eval (from likelihood) bdk[%d]terms(%lf) = %lf", t,b,val);
#endif
  return val;
}
#else
static double bdkterms(double b, void *mydata) {
  int t;
  double val;
  for (t=0; t<ddN.T; t++) 
    ddP.bdk[t] = b;
  val = likelihood_bdk_asym(-1);
  myarms_evals++;
#ifdef B_DEBUG
  yap_message("Eval (from likelihood) bdkterms(%lf) = %lf\n", b, val);
#endif
  return val;
}
#endif

/*
 *   same logic as likelihood_bdk() but we store stats
 */
static void bdk_store() {
  D_MiSi_t dD;   
  int t, i;
  int totmemory = 0;

  misi_init(&ddM,&dD);
  docstats = malloc(sizeof(*docstats)*ddN.DT);
  totmemory += sizeof(*docstats)*ddN.DT;
  if ( !docstats )
    yap_quit("Cannot allocate memory in bdk_store()\n");

  /*
   *    do for each document in turn;
   *    have to build the counts and table counts for
   *    each from the t,r values in ddS.z[]
   */
  for (i=0; i<ddN.DT; i++) {
    int docmemory;
    misi_build(&dD, i, 0);
    /*
     *  get memory usage for this doc
     */
    docmemory = 1;
    for (t=0; t<ddN.T; t++) {
      if ( dD.Mi[t]>0 ) {
	docmemory += 3;
      }
    }
    if ( docmemory==1 ) {
      docstats[i] = NULL;
      continue;
    }
    docstats[i] = u16vec(docmemory);
    totmemory += docmemory*sizeof(docstats[i][0]);
    {
      int mem = 0;
      uint16_t *dstats = docstats[i];
      for (t=0; t<ddN.T; t++) {
	if ( dD.Mi[t]>0 ) {
	  dstats[mem] = t;
	  dstats[mem+1] = dD.Mi[t];
	  dstats[mem+2] = dD.Si[t];
	  mem += 3;
	}
      }
      dstats[mem] = ddN.T+1;
      assert(mem+1==docmemory);
    }
    misi_unbuild(&dD, i, 0);
  }
  misi_free(&dD);
  if ( verbose>1 )
    yap_message("bdk_store:  built BDK stats in %d bytes\n", totmemory);
}
static void bdk_free() {
  int i;
  for (i=0; i<ddN.DT; i++) 
    if ( docstats[i]!=NULL )
      free(docstats[i]);
  free(docstats);
}


/*
 *    this uses ARMS
 */
void sample_bdk(double *b) {
#ifdef SAMPLEBDK_MULTI
  static int batchstart = 0;
  int tt;
#endif
  double startlike;
  int t;
  assert(b);
  bdk_store();
#ifdef SAMPLEBDK_MULTI
  assert(ddP.bdkbatch>0);
  if ( verbose>1 ) {
      t = batchstart%ddN.T;
      startlike = likelihood();
      yap_message("sample_bdk[%d:%d] (pre): b=%lf, lp=%lf\n",
		  t, (ddP.bdkbatch+batchstart-1)%ddN.T,b[t], startlike);
  }
  for (tt=0; tt<ddP.bdkbatch; tt++) {
    t = (batchstart+tt)%ddN.T;
    
    myarmsMH(PYP_CONC_MIN, PYP_CONC_MAX, &bdkterms, &t, &b[t], "bdk", 1);
  }
  if ( verbose>1 ) {
      double endlike = likelihood();
      t = (batchstart)%ddN.T;
      yap_message("sample_bdk[%d:%d] (post): bdk=%lf, lp=%lf\n",
		  t, (ddP.bdkbatch+batchstart-1)%ddN.T, b[t], endlike);
      if ( endlike < startlike-50 ) {
	yap_quit("Sampler failed due to huge decrease!\n");
      }
  }
  batchstart = (batchstart+tt)%ddN.T;
#else    
  if ( verbose>1 ) {
    startlike = likelihood();
    yap_message("sample_bdk (pre): bdk[*]=%lf, lp=%lf\n",
		b[0], startlike);
  }
    
  myarmsMH(PYP_CONC_MIN, PYP_CONC_MAX, &bdkterms, NULL, &b[0], "bdk", 1);
  for (t=1; t<ddN.T; t++)
    b[t] = b[0];
  if ( verbose>1 ) {
    double endlike = likelihood();
    yap_message("sample_bdk (post): bdk[*]=%lf, lp=%lf\n",
		b[0], endlike);
    if ( endlike < startlike-50 ) {
      yap_quit("Sampler failed due to huge decrease!\n");
    }
  }
#endif
  bdk_free();
}

 
