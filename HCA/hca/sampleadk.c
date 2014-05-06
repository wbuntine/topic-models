/*
 * Sampling utility for adk
 * Copyright (C) 2011-2013 Wray Buntine 
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
#include "stable.h"
#include "sample.h"
#include "stats.h"

//  #define A_DEBUG

#ifdef A_DEBUG
  static double last_val = 0;
  static double last_like = 0;
#endif
/*
 *  globals from hca.h and stats.h
 */
extern int verbose;
double likelihood_bdk();
  
/*
 *    docstats[d] points to store for doc d
 *    
 *    stores blocks per nonzero topic
 *       [0] = topic
 *       [1] = total count
 *       [2] = total tables
 *       [3] = no. of Mi[] Si[] pairs following with Mi[]>=2
 *       [4+2w]+[5+2w] = Mi[]+Si[] for w-th word
 *
 *    ending block
 *       [0] = ddN.T+1
 */
static uint16_t **docstats;

/*
 *     uset<0     ==>   compute for all t
 *     otherwise  ==>   compute for just that t
 */
static double likelihood_ad_asym() {
  double *lgabd = dvec(ddN.T);
  double la;
  int t, i;
  double likelihood = 0;
  la = log(ddP.ad);
  for (t=0; t<ddN.T; t++) {
    lgabd[t] = lgamma(ddP.bdk[t]/ddP.ad);
  }
  for (i=0; i<ddN.DT; i++) {
    uint16_t *store = docstats[i];
    int mem3, j;
    if ( store==NULL ) 
	continue;
    while ( store[0]<ddN.T ) {
      t = store[0];
      assert(store[1]>0);
      likelihood += ((int)store[2])*la 
	+ gammadiff((int)store[2], ddP.bdk[t]/ddP.ad, lgabd[t]);
      mem3 = store[3];
      store += 4;
      for (j=0; j<mem3; j++) {
	likelihood += S_S(ddC.SD,store[j*2],store[j*2+1]);
      }
      store += 2*mem3;
    }
    yap_infinite(likelihood);
  }
  free(lgabd);
  return likelihood;
}


/*
 */
static double adkterms(double mya, void *mydata) {
  double val = 0;
#ifdef A_DEBUG
  float save_a = ddC.SD->a;
  double like;
#endif
  ddP.ad = mya;
  cache_update("ad");
  val = likelihood_ad_asym();
  myarms_evals++;
#ifdef A_DEBUG
  yap_message("Eval adkterms(%lf) = %lf", mya, val);
  like = likelihood_bdk();
  if ( last_val != 0 ) {
    yap_message(", lp=%lf diffs=%lf vs %lf\n", like, 
		val-last_val, like-last_like);
  } else
    yap_message("\n");
  last_like = like;
  last_val = val;
#endif
  return val;
}


/*
 *   same logic as likelihood_bdk() but we store stats
 */
static void adk_store() {
  D_MiSi_t dD;
  int mi, l, t, i;
#define MAXW 10
  /*
   *    store[t][0] = #words with Mi[word][t]>1
   *    store[t][1+2w] = Mi[word][t]
   *    store[t][2+2w] = Si[word][t]
   */
  uint16_t **store;
  uint16_t **overflow;
  int totmemory = 0;

  misi_init(&ddM,&dD);
  store = u16mat(ddN.T,MAXW*2+1);
  overflow = malloc(ddN.T*sizeof(overflow));
  docstats = malloc(sizeof(*docstats)*ddN.DT);
  totmemory += sizeof(*docstats)*ddN.DT;
  if ( !docstats  || !overflow)
    yap_quit("Cannot allocate memory in bdk_store()\n");
  for (t=0; t<ddN.T; t++) 
    overflow[t] = NULL;

  /*
   *    do for each document in turn;
   *    have to build the counts and table counts for
   *    each from the t,r values in ddS.z[]
   */

  for (i=0; i<ddN.DT; i++) {
    int docmemory;
    misi_build(&dD, i, 0);

    /*  compute likelihood and zero too*/
    mi = ddM.MI[i];      
    for (l=ddD.NdTcum[i]; l< ddD.NdTcum[i+1]; l++) {
      if ( M_multi(l) ) {
	int mii = ddM.multiind[mi] - dD.mi_base;
	t = Z_t(ddS.z[l]);
	if ( dD.Mik[mii][t]>0 ) {
	  if ( dD.Mik[mii][t]>1 ) {
	    int scnt = store[t][0];
	    /*
	     *    save greater than 2 counts in store+overflow
	     */
	    if ( scnt>=MAXW ) {
	      if ( scnt==MAXW ) {
		assert(overflow[t]==NULL);
		overflow[t] = malloc(2*MAXW*sizeof(overflow[t][0]));
	      } else if ( scnt%MAXW == 0 ) {
		int mysize = ((scnt-MAXW)/MAXW) + 1;
		overflow[t] = realloc(overflow[t],
				      2*mysize*MAXW*sizeof(overflow[t][0]));
	      }
	      if ( !overflow[t] )
		yap_quit("Cannot allocate memory in bdk_store()\n");
	      overflow[t][2*(scnt-MAXW)] = dD.Mik[mii][t];
	      overflow[t][1+2*(scnt-MAXW)] = dD.Sik[mii][t];
	    } else {
	      store[t][1+2*scnt] = dD.Mik[mii][t];
	      store[t][2+2*scnt] = dD.Sik[mii][t];
	    }
	    store[t][0]++;
	  }
	  /*  zero these now so don't double count later */
	  dD.Mik[mii][t] = 0;
	  dD.Sik[mii][t] = 0;
	}
	mi++;
      }
    }

    /*
     *  get memory usage for this doc
     */
    docmemory = 1;
    for (t=0; t<ddN.T; t++) {
      if ( dD.Mi[t]>0 ) {
        docmemory += 4 + store[t][0]*2;
      }
    }
    if ( docmemory==1 ) {
      docstats[i] = NULL;
      continue;
    }
    docstats[i] = u16vec(docmemory);
    totmemory += docmemory*sizeof(docstats[i][0]);
    {
      int mem = 0, j;
      uint16_t *dstats = docstats[i];
      for (t=0; t<ddN.T; t++) {
	if ( dD.Mi[t]>0 ) {
	  dstats[mem] = t;
	  dstats[mem+1] = dD.Mi[t];
	  dstats[mem+2] = dD.Si[t];
	  dstats[mem+3] = store[t][0];
	  mem += 4;
	  for (j=0; j<store[t][0] && j<MAXW; j++) {
	    dstats[mem+2*j] = store[t][1+2*j];
	    dstats[mem+2*j+1] = store[t][2+2*j];
	    store[t][1+2*j] = 0;
	    store[t][2+2*j] = 0;
	  }
	  for ( ; j<store[t][0]; j++) {
	    dstats[mem+2*j] = overflow[t][2*(j-MAXW)];
	    dstats[mem+2*j+1] = overflow[t][1+2*(j-MAXW)];
	  }
	  mem += store[t][0]*2;
	  if ( overflow[t] )
	    free(overflow[t]);
	  overflow[t] = NULL;
	  store[t][0] = 0;
        }
      }
      dstats[mem] = ddN.T+1;
      assert(mem+1==docmemory);
    }
    misi_unbuild(&dD, i, 0);
  }
  misi_free(&dD);
  free(overflow);
  free(store[0]);  free(store);
  if ( verbose>1 )
    yap_message("bdk_store:  built BDK stats in %d bytes\n", totmemory);
}

static void adk_free() {
  int i;
  for (i=0; i<ddN.DT; i++)
    if ( docstats[i]!=NULL )
      free(docstats[i]);
  free(docstats);
}


void sample_adk(double *mya) {
#ifdef A_DEBUG
  last_val = 0;
  last_like = 0;
#endif
  adk_store();
  if ( verbose>1 )
    yap_message("sample_adk (pre):  ad=%lf, lp=%lf\n",
		*mya, likelihood());
  myarms(PYP_DISC_MIN, PYP_DISC_MAX, &adkterms, NULL, mya, "adk");
  cache_update("ad");
  adk_free();
  if ( verbose>1 )
    yap_message("sample_adk (post):  ad=%lf, lp=%lf\n",
		*mya, likelihood());
 }

