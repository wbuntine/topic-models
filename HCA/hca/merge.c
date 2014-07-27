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
#include "heap.h"
#include "stats.h"
#include "cache.h"

/*
 *    adapted by diffing code in like.c
 */ 

static double likemerge_DIRalpha(k1,k2) {
  /*
   *   Dirichlet for topics
   */
  int i;
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
  int j;
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

int fveccmp(uint32_t k1, uint32_t k2, void *par) {
  float *fvec = (float *)par;
  if ( fvec[k1]>=fvec[k2] ) 
    return 1;
  return 0;
}

/*
 *  place all topic k2 data in topic k1;
 *  find optimal table counts Tdt[.][k1];
 *  return likelihood ratio for this;
 *  Tdt[] must be assigned on input
 *        Tdt = u16vec(ddN.DT); 
 */
double likemerge_alpha(int k1, int k2, uint16_t *Tdt) {
  int d, k;

  /*
   *   storing/saving the new/proposed version
   */
  uint32_t TDt, TDTm;
  /*
   *  cache for doc side change
   */
  uint16_t *TdT;
  uint16_t *Ndt;

  struct heap_s up;
  struct heap_s down;
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

  double likelihood;

  if ( k1<=0 || k2<=0 || ddS.TDt[k1]==0 || ddS.TDt[k2]==0 ) 
    return 0.0;
 
  /*
   *   build local store
   */
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
  TDTm = ddS.TDT - TDt;

  /*
   *  initialise sort
   */
  for (d=0; d<ddN.DT; d++) {
    /*   don't change for some docs */
    Tdt_up[d] = d;
    Tdt_down[d] = d;
    if ( Tdt[d]<Ndt[d] )
      score_up[d] = (ddP.bpar + ddP.apar*TdT[d]) / S_V(ddC.SX,Ndt[d],Tdt[d]+1);
    else 
      score_up[d] = 0;
    if ( Tdt[d]>1 )
      score_down[d] = S_V(ddC.SX,Ndt[d],Tdt[d])/(ddP.bpar + ddP.apar*(TdT[d]-1));
    else
      score_down[d] = 0;
  }
  assert(TDt>0);

  /*  
   *  use a heap, so only top of heap is least 
   */
  heap_init(&up, Tdt_up, ddN.DT, fveccmp, (void *)score_up);
  heap_init(&down, Tdt_down, ddN.DT, fveccmp, (void *)score_down);

  while ( 1 ) {
    float upv;
    float downv;
    upv = merge_alphabasetopicprob(TDTm+TDt, TDt, k1)*score_up[heap_front(&up)];
    if ( TDt>1 )
      downv = score_up[heap_front(&up)] 
        / merge_alphabasetopicprob(TDTm+TDt-1, TDt-1, k1);
    else
      downv = 0.0;
    if ( downv>upv && downv>1.0 ){
      //  decrement this
      d = heap_front(&down);
      TdT[d]--;
      Tdt[d]--;
      assert(Tdt[d]>0);
      TDt--;
      heap_pop(&down);
      heap_remove(&up,d);
    } else if ( downv<upv && upv>1.0 ){
      //  increment this
      d = heap_front(&up);
      TdT[d]++;
      Tdt[d]++;
      assert(Tdt[d]<=Ndt[d]);
      TDt++;
      heap_pop(&up);
      heap_remove(&down,d);
    } else {
      //  none are better so quit
      break;
    }
    if ( Tdt[d]<Ndt[d] )
      score_up[d] = (ddP.bpar + ddP.apar*TdT[d]) / S_V(ddC.SX,Ndt[d],Tdt[d]+1);
    else 
      score_up[d] = 0;
    if ( Tdt[d]>1 )
      score_down[d] = S_V(ddC.SX,Ndt[d],Tdt[d])/(ddP.bpar + ddP.apar*(TdT[d]-1));
    else
      score_down[d] = 0;
    /*
     *  now adjust the two heaps for new vals for [d]
     */
    heap_push(&down,d);
    heap_push(&up,d);
  }
  
  /*
   *  so have optimal Tdt[] for merge;
   *  compute final likelihood ratio based on this
   */
  {
    double la = 0;
    double lb = log(ddP.bpar);
    int TD_diff = 0;
    if ( ddP.apar>0 ) la = log(ddP.apar);
    likelihood = 0;
    for (d=0; d<ddN.DT; d++) {
      int Td_diff;
      if ( Ndt[d]>1 ) {
        likelihood -= S_S(ddC.SX,ddS.Ndt[d][k2],ddS.Tdt[d][k2]);
        likelihood -= S_S(ddC.SX,ddS.Ndt[d][k1],ddS.Tdt[d][k1]);
        likelihood += S_S(ddC.SX,Ndt[d],Tdt[d]);
      }
      yap_infinite(likelihood);
      TD_diff += Td_diff = (Tdt[d]-ddS.Tdt[d][k2]-ddS.Tdt[d][k1]);
      if ( Td_diff==0 )
        continue;
      if ( ddP.apar==0 ) {
        likelihood += Td_diff*lb;
      } else {
        likelihood += Td_diff*la;
        if ( Td_diff<0 ) 
          likelihood -= gammadiff(-Td_diff, TDt-1+ddP.bpar/ddP.apar, 0.0);
        else
          likelihood += gammadiff(Td_diff, TDt-Td_diff-1+ddP.bpar/ddP.apar, 0.0);
      }
      yap_infinite(likelihood);
    }      
    if ( TD_diff!=0 ) {
      if ( ddP.PYalpha==H_PDP ) {
        likelihood += TD_diff*log(ddP.alphapr[k1]);
      } else if ( ddP.PYalpha==H_HDP ) {
        likelihood += lgamma(TDTm+TDt-TD_diff+ddP.b0) - lgamma(TDTm+TDt+ddP.b0);
        if ( ddS.TDt[k2]>0 )
          likelihood -= gammadiff(ddS.TDt[k2], ddP.b0*ddP.alphapr[k2], 0.0);
        if ( TD_diff )
          likelihood += gammadiff(TD_diff, TDt-TD_diff+ddP.b0*ddP.alphapr[k1], 0.0);
        else
          likelihood -= gammadiff(-TD_diff, TDt+ddP.b0*ddP.alphapr[k1], 0.0);
      } else {
        double lga0 = lgamma(1-ddP.a0);
        likelihood += lgamma(TDTm+TDt-TD_diff+ddP.b0) - lgamma(TDTm+TDt+ddP.b0);
        if ( ddS.TDt[k2]>0 )
          likelihood -= log(ddP.b0+ddP.a0*(ddS.TDTnz-1));
        if ( ddS.TDt[k2]>1 )
          likelihood -= lgamma(ddS.TDt[k2]-ddP.a0) - lga0;
        if ( ddS.TDt[k1]>1 )
          likelihood -= lgamma(ddS.TDt[k1]-ddP.a0) - lga0;
        likelihood += lgamma(TDt-ddP.a0) - lga0;
      }
    }
  }

  free(TdT);
  free(Ndt);
  free(Tdt);
  free(score_up);
  free(score_down);
  heap_free(&up);
  heap_free(&down);
  yap_infinite(likelihood);
  return likelihood;
}


