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

double likemerge_DIRalpha(k1,k2) {
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

double likemerge_DIRbeta(k1,k2) {
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

double alphabasechange(int t) {
  if ( ddP.PYalpha==H_HDP )
    return (TDt + ddP.b0*ddP.alphapr[t]);
  if ( ddP.PYalpha==H_PDP )
    return ddP.alphapr[t];
  return (TDt - ddP.a0)/(TDT + ddP.b0)
}

/*
 *  place all topic k2 data in topic k1
 */
double likemerge(int k1, int k2) {
  int t, d, k;
  uint16_t *Tdt;
  uint16_t *TdT;
  uint16_t *Tdt_up;
  uint16_t *Tdt_down;
  uint32_t TDt;
  float *score_up;
  float *score_down;
  float base_up;
  float base_down;
  double likelihood = 0;

  if ( k1<=0 || k2<=0 ) 
    return 0.0;
 
  TdT = u16vec(ddN.DT);
  Ndt = u16vec(ddN.DT);
  Tdt = u16vec(ddN.DT);
  Tdt_up = u16vec(ddN.DT);
  Tdt_down = u16vec(ddN.DT);
  score_up = u16vec(ddN.DT);
  score_down = u16vec(ddN.DT);
  if ( !score_down || !score_up || !Ndt || !TdT ||
       !Tdt || !Tdt_up || !Tdt_down )
    yap_quit("Out of memory in likemerge()\n");

  TDt = 0;
  for (d=0; d<ddN.DT; d++) {
    Ndt[d] = ddS.Ndt[d][k1] + ddS.Ndt[d][k1];
    TdT[d] = 0;  
    for (k=0; k<ddN.T; k++)
      TdT[d] += ddS.Tdt[d][k];
    Tdt[d] = ddS.Tdt[d][k1] + ddS.Tdt[d][k1];
    if ( Tdt[d]>2 )
      Tdt[d]--;
    TDt += Tdt[d];
  }
  for (d=0; d<ddN.DT; d++) {
    Tdt_up[d] = d;
    Tdt_down[d] = d;
    score_up[d] = (ddP.bpar + ddP.apar*TdT[d])
      * S_V(ddP.SX,Ndt[d],Tdt[d]+1) * alphabasetopicprob*(t);
    score_down[d] = 1.0/(ddP.bpar + ddP.apar*(TdT[d]-1))
      / S_V(ddP.SX,Ndt[d],Tdt[d]) / alphabasetopicprob*(t-1);
    ????;
  }

  free(score_up);
  free(score_down);
  free(TdT);
  free(Tdt);
  free(Tdt_up);
  free(Tdt_down);
  yap_infinite(likelihood);
  return likelihood;
}


