/*
 * Likelihood diff due to merge
 * Copyright (C) 2014 Wray Buntine
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@monash.edu)
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "yap.h"
#include "util.h"
#include "stable.h"
#include "lgamma.h"
#include "hca.h"
#include "stats.h"
#include "cache.h"

/*
 *    adapted by diffing code in like.c
 */ 

static double likemerge_DIRalpha(k1,k2) {
  /*
   *   Dirichlet for topics
   */
  int i,t;
  double likelihood = 0;
  for (i=0; i<ddN.DT; i++) {
    if ( ddS.Ndt[i][k2]>0 ) {
      likelihood +=
	gammadiff((int)ddS.Ndt[i][k1]+ddS.Ndt[i][k2], ddP.alphapr[k1], 0.0)
	- gammadiff((int)ddS.Ndt[i][k1], ddP.alphapr[k1], 0.0)
	- gammadiff((int)ddS.Ndt[i][k2], ddP.alphapr[k2], 0.0);
    }
  }
  yap_infinite(likelihood);
  return likelihood;
}

static double likemerge_DIRbeta(k1,k2) {
  int j,t;
  double likelihood = 0;
  double val = 0;
  for (j=0; j<ddN.W; j++) {
    if ( ddS.Nwt[j][k2]>0 ) {
      assert(ddP.betapr[j]>0);
      val += gammadiff((int)ddS.Nwt[j][k1]+ddS.Nwt[j][k2], ddP.betapr[j], 0.0)
	- gammadiff((int)ddS.Nwt[j][k1], ddP.betapr[j], 0.0)
	- gammadiff((int)ddS.Nwt[j][k2], ddP.betapr[j], 0.0);
    } 
  }     
  val -= gammadiff((int)ddS.NWt[k1]+ddS.NWt[k2], ddP.betatot, 0.0)
    - gammadiff((int)ddS.NWt[k1], ddP.betatot, 0.0)
    - gammadiff((int)ddS.NWt[k2], ddP.betatot, 0.0);
  likelihood += val;
  yap_infinite(likelihood);
  return likelihood;
}

static double merge_alphabasetopicprob(uint32_t TDT, uint32_t TDt, int t) {
  assert(t>=0);
  if ( ddP.PYalpha==H_PDP ) 
    return ddP.alphapr[t];
  if ( ddP.PYalpha==H_HDP ) {
    assert(ddP.a0==0);
    return  ( (double)TDt+ddP.b0*ddP.alphapr[t] )
      /((double)TDT+ddP.b0);
  }
  /*
   *   normal situation when empty topics exist
   */
  return  ((double)TDt-ddP.a0)/((double)TDT+ddP.b0);	 
}

/*
 *  place all topic k2 data in topic k1
 */
double likemerge(int k1, int k2) {
  int t, d, k;

  /*
   *   storing/saving the new/proposed version
   */
  uint32_t TDt;
  uint16_t *Tdt;
  /*
   *  cache for doc side change
   */
  uint16_t *TdT;
  uint16_t *Ndt;

  int sortD;
  /*
   *    sorting on moves,
   *    stores doc index for 
   */
  uint32_t *Tdt_up;
  uint32_t *Tdt_down;
  /*
   *     change from incr/decr this docs Tdt
   */
  float *score_up;
  float *score_down;
  /*
   *     change from incr/decr total TDt
   */
  float base_up;
  float base_down;

  double likelihood = 0;

  if ( k1<=0 || k2<=0 || ddS.TDt[k1]==0 || ddS.TDt[k2]==0 ) 
    return 0.0;
 
  /*
   *   build local store
   */
  Tdt = u16vec(ddN.DT); 
  TdT = u16vec(ddN.DT); 
  Ndt = u16vec(ddN.DT);
  Tdt_up = u32vec(ddN.DT);
  Tdt_down = u32vec(ddN.DT);
  score_up = fvec(ddN.DT);
  score_down = fvec(ddN.DT);
  if ( !score_down || !score_up || !Ndt || !TdT ||
       !Tdt || !Tdt_up || !Tdt_down )
    yap_quit("Out of memory in likemerge()\n");

  /*
   *  initialise all data entries to simple sum
   */
  TDt = 0;
  for (d=0; d<ddN.DT; d++) {
    Tdt[d] = ddS.Tdt[d][k1] + ddS.Tdt[d][k2];
    Ndt[d] = ddS.Ndt[d][k1] + ddS.Ndt[d][k2];
    TdT[d] = 0;
    for (k=0; k<ddN.T; k++)
      TdT[d] += ddS.Tdt[d][k];
    TDt += Tdt[d];
  }

  /*
   *  initialise sort
   */
  sortD = 0;
  for (d=0; d<ddN.DT; d++) {
    /*   don't change for some docs */
    if ( Ndt[d]<=1 )
      continue;
    Tdt_up[sortD] = d;
    Tdt_down[sortD] = d;
    sortD++;
    if ( Tdt[d]<Ndt[d] )
      score_up[d] = (ddP.bpar + ddP.apar*TdT[d]) / S_V(ddP.SX,Ndt[d],Tdt[d]+1) 
	* merge_alphabasetopicprob(TDT, TDt, t);
    else 
      score_up[d] = 0;
    if ( Tdt[d]>1 )
      score_down[d] = S_V(ddP.SX,Ndt[d],Tdt[d])/(ddP.bpar + ddP.apar*(TdT[d]-1))
	/ merge_alphabasetopicprob(TDT-1, TDt-1, t);
    else
      score_down[d] = 0;
  }
  /*  sortD is count of docs to be optimised */
  /*  use a heap, so only top of heap is least */

  free(TdT);
  free(Ndt);
  free(Tdt);
  free(score_up);
  free(score_down);
  free(Tdt_up);
  free(Tdt_down);
  yap_infinite(likelihood);
  return likelihood;
}


