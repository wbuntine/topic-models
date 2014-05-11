/*
 * Deal with the caches and tables
 * Copyright (C) 2012-2014 Wray Buntine 
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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "yap.h"
#include "yaps.h"
#include "lgamma.h"
#include "hca.h"
#include "data.h"
#include "pctl.h"
#include "cache.h"

#define mymax(A,B) ((A>B)?(A):(B))

#ifdef CACHE_ABTP
double alphabasetopicprob(int t);
#endif

void cache_init(int maxM, int maxW) {
  /*
   *   now the libstb library has its own separate error handling,
   *   so make sure it goes through our yapper
   */
  yaps_yapper(yap_va);

  /*
   *   we'll share tables for all those with a==0,
   *   since it never changes
   */
  {
    char tag[200];
    int mN=0, mM=0, sN=0, sM=0;
    ddC.S0 = NULL;
    tag[0] = 0;
    if ( ddP.apar==0 ) {
      int maxD = ddD.NdTmax*1.5;
      mM = mymax(mM,maxM);
      mN = mymax(mN,maxD+1);
      sM = mymax(sM,100);
      sN = mymax(sN,maxD+1);
      strcat(tag, "a, ");
    }
    if ( ddP.ad==0 ) {
      int nmax=ddM.Mi_max+10;
      mM = mymax(mM,nmax);
      mN = mymax(mN,nmax);
      sM = mymax(sM,nmax);
      sN = mymax(sN,nmax);
      strcat(tag, "ad, ");
    }
    if ( ddP.awpar==0 ) {
      mM = mymax(mM,maxM);
      mN = mymax(mN,maxW);
      sM = mymax(sM,1000);
      sN = mymax(sN,100);
      strcat(tag, "aw, ");
    }
    if ( mM>0 ) {
      ddC.S0 = S_make(sN, sM, mN, mM, 0,
#ifdef H_THREADS
                      S_THREADS|
#endif
                      S_STABLE|S_UVTABLE|S_FLOAT| S_QUITONBOUND);
      if ( ddP.apar==0 ) 
        ddC.SX = ddC.S0;
      if ( ddP.awpar==0 )
        ddC.SY = ddC.S0;
      if ( PCTL_BURSTY() && ddP.ad==0 ) 
        ddC.SD = ddC.S0;
      strcat(tag, " all zero PYP");
      S_tag(ddC.S0, tag);
      if ( verbose ) 
        S_report(ddC.S0,NULL);
    }
  }
  if ( ddP.bdk!=NULL && ddP.ad!=0 ) {
    /*
     *   clear maximum is number times any given word appears in
     *   any one doc (max for any word and any doc)
     */
    int nmax=ddM.Mi_max+10;
    ddC.SD = S_make(nmax, nmax, nmax, nmax, ddP.ad, 
#ifdef H_THREADS
 		    S_THREADS|
#endif
		    S_STABLE|S_UVTABLE|S_FLOAT| S_QUITONBOUND);
    S_tag(ddC.SD,"SD, doc PYP");
    if ( verbose ) S_report(ddC.SX,NULL);
  }
  if ( ddP.PYalpha ) {
    /*
     *    clear maximum of N is maximum number of words;
     *    in a document use the passed in limit for T
     */
    int maxD = ddD.NdTmax*1.5;
#ifdef CACHE_ABTP
    alphabasetopicprob(-1);
#endif
    if ( ddP.apar!=0 ) {
      ddC.SX = S_make(maxD+1, 100, maxD+1, maxM, ddP.apar, 
#ifdef H_THREADS
                      S_THREADS|
#endif
                      S_STABLE|S_UVTABLE|S_FLOAT| S_QUITONBOUND);
      S_tag(ddC.SX,"SX, docXtopic PYP");
      if ( verbose ) S_report(ddC.SX,NULL);
    }
    gcache_init(&ddC.qda, ddP.apar);
    gcache_init(&ddC.lgb, ddP.bpar);
    if ( ddP.apar>0 )
      gcache_init(&ddC.lgba, ddP.bpar/ddP.apar);
  } else {
    gcache_init(&ddC.lgalpha, ddP.alpha);
    gcache_init(&ddC.lgtotalpha, ddN.T*ddP.alpha);
  }
  if ( ddP.PYbeta ) {
    if ( ddP.awpar!=0 ) {
      ddC.SY = S_make(1000, 100, maxW, maxM, ddP.awpar,
#ifdef H_THREADS
                      S_THREADS|
#endif
                      S_STABLE|S_UVTABLE|S_FLOAT| S_QUITONBOUND);
      S_tag(ddC.SY,"SY, topicXword PYP");
      if ( verbose ) S_report(ddC.SY,NULL);
    }
    gcache_init(&ddC.qdaw, ddP.awpar);
    gcache_init(&ddC.lgbw, ddP.bwpar);
    if ( ddP.awpar>0 )
      gcache_init(&ddC.lgbaw, ddP.bwpar/ddP.awpar);
  } else {
    if ( ddP.betac>0 )
      gcache_init(&ddC.lgbetac, ddP.betac);
    gcache_init(&ddC.lgbeta, ddP.beta);
  }
}

void cache_free() {
#ifdef CACHE_ABTP
  alphabasetopicprob(-5*ddN.T);
#endif
  if ( ddP.PYbeta &&  ddP.awpar!=0 ) {
    S_free(ddC.SY);
  }  
  if ( ddP.PYalpha &&  ddP.apar!=0 ) {
    S_free(ddC.SX);
  }
  if ( ddP.bdk!=NULL && ddP.ad!=0 ) {
    S_free(ddC.SD);
  }
}

void cache_update(char *par) {
  if ( ddP.bdk!=NULL ) {
    if ( strcmp(par,"ad")==0 ) {
      if ( ddC.SD->a != ddP.ad )
	S_remake(ddC.SD,ddP.ad);
    } 
  }
  if ( ddP.PYalpha ) {
    if ( strcmp(par,"b")==0 ) {
      gcache_init(&ddC.lgb, ddP.bpar);
      if ( ddP.apar>0 )
	gcache_init(&ddC.lgba, ddP.bpar/ddP.apar);
    } else if ( strcmp(par,"a")==0 ) {
      gcache_init(&ddC.lgba, ddP.bpar/ddP.apar);
      gcache_init(&ddC.qda, ddP.apar);
      if ( ddC.SX->a != ddP.apar )
	S_remake(ddC.SX,ddP.apar);
    } 
#ifdef CACHE_ABTP
    else if ( strcmp(par,"b0")==0 || strcmp(par,"a0")==0 ) {
      alphabasetopicprob(-(ddN.T+1));
    }
#endif
  } else if ( strcmp(par,"alpha")==0 ) { 
    gcache_init(&ddC.lgalpha, ddP.alpha);
    gcache_init(&ddC.lgtotalpha, ddN.T*ddP.alpha);
  }
  if ( ddP.PYbeta ) {
    if ( strcmp(par,"bw")==0 ) { 
      gcache_init(&ddC.lgbw, ddP.bwpar);
      if ( ddP.awpar>0 )
	gcache_init(&ddC.lgbaw, ddP.bwpar/ddP.awpar);
    } else if ( strcmp(par,"aw")==0 ) { 
      gcache_init(&ddC.lgbaw, ddP.bwpar/ddP.awpar);
      gcache_init(&ddC.qdaw, ddP.awpar);
      if ( ddC.SY->a != ddP.awpar ) 
	S_remake(ddC.SY, ddP.awpar);
    } 
  } else {
    if ( strcmp(par,"beta")==0 ) {
      fixbeta(NULL, NULL);
      if ( ddP.betac>0) {
	gcache_init(&ddC.lgbeta, ddP.betac);
      }
      gcache_init(&ddC.lgbeta, ddP.beta);
    }
  }
}
