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

//   switch off TXt[] optimisation
#define NOOPT_MERGE

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

#ifndef NOOPT_MERGE
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

static int fveccmp(uint32_t k1, uint32_t k2, void *par) {
  float *fvec = (float *)par;
  if ( fvec[k1]>=fvec[k2] ) 
    return 1;
  return 0;
}
#endif
/*
 *  place all topic k2 data in topic k1;
 *  find optimal table counts Tdt[.][k1];
 *  return likelihood ratio for this;
 *  Tdt[] must be assigned on input
 *        Tdt = u16vec(ddN.DT); 
 */
typedef struct merge_s {
  /*
   *   storing/saving the new/proposed version
   */
  uint32_t TDt, TDTm;
  /*
   *  cache for doc side change
   */
  uint16_t *TdT;
  uint16_t *Ndt;
  uint16_t *Tdt;
} merge_t;


static void merge_free(merge_t *M) {
  free(M->Ndt);
  free(M->TdT);
  free(M->Tdt);
}

static void merge_init(int k1, int k2, merge_t *M) {
  int d, k;
  /*
   *   build local store
   */
  M->Tdt = u16vec(ddN.DT); 
  M->TdT = u16vec(ddN.DT); 
  M->Ndt = u16vec(ddN.DT);
  if ( !M->Ndt || !M->TdT || !M->Tdt  )
    yap_quit("Out of memory in likemerge()\n");

  /*
   *  initialise all data entries to simple sum
   */
  M->TDt = 0;
  for (d=0; d<ddN.DT; d++) {
    M->Tdt[d] = ddS.Tdt[d][k1] + ddS.Tdt[d][k2];
    M->Ndt[d] = ddS.Ndt[d][k1] + ddS.Ndt[d][k2];
    M->TdT[d] = 0;
    for (k=0; k<ddN.T; k++)
      M->TdT[d] += ddS.Tdt[d][k];
    M->TDt += M->Tdt[d];
    assert(M->Tdt[d]<=M->Ndt[d]);
    if ( d>0 )
      assert(M->Tdt[d-1]<=M->Ndt[d-1]);
  }
  M->TDTm = ddS.TDT - M->TDt;
#ifndef NDEBUG
  for (d=0; d<ddN.DT; d++) {
    assert(M->Tdt[d]<=M->Ndt[d]);
  }
#endif
}

#ifndef NOOPT_MERGE
static void merge_opt_Tdt(int k1, int k2, merge_t *M) {
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

  Tdt_up = u32vec(ddN.DT);
  Tdt_down = u32vec(ddN.DT);
  score_up = fvec(ddN.DT);
  score_down = fvec(ddN.DT);
  if ( !score_down || !score_up || !Tdt_up || !Tdt_down )
    yap_quit("Out of memory in likemerge()\n");
  /*
   *  initialise sort
   */
  for (d=0; d<ddN.DT; d++) {
    assert(M->Tdt[d]<=M->Ndt[d]);
    /*   don't change for some docs */
    Tdt_up[d] = d;
    Tdt_down[d] = d;
    if ( M->Tdt[d]<M->Ndt[d] )
      score_up[d] = (ddP.bpar + ddP.apar*M->TdT[d]) 
	/ S_V(ddC.SX,M->Ndt[d],M->Tdt[d]+1);
    else 
      score_up[d] = 0;
    if ( M->Tdt[d]>1 )
      score_down[d] = S_V(ddC.SX,M->Ndt[d],M->Tdt[d])
	/(ddP.bpar + ddP.apar*(M->TdT[d]-1));
    else
      score_down[d] = 0;    
    assert((M->Tdt[d]>1)||score_down[d]==0);
    assert((M->Tdt[d]<M->Ndt[d])||score_up[d]==0);
    assert(M->Tdt[d]<=M->Ndt[d]);
  }
  assert(M->TDt>0);

  /*  
   *  use a heap, so only top of heap is least 
   */
  heap_init(&up, Tdt_up, ddN.DT, fveccmp, (void *)score_up);
  heap_init(&down, Tdt_down, ddN.DT, fveccmp, (void *)score_down);

  while ( 1 ) {
    float upv;
    float downv;
    upv = merge_alphabasetopicprob(M->TDTm+M->TDt, M->TDt, k1)
      *score_up[heap_front(&up)];
    if ( M->TDt>1 )
      downv = score_down[heap_front(&down)] 
        / merge_alphabasetopicprob(M->TDTm+M->TDt-1, M->TDt-1, k1);
    else
      downv = 0.0;
    if ( downv>upv && downv>1.0 ){
      //  decrement this
      d = heap_front(&down); 
      M->TdT[d]--;
      M->Tdt[d]--;
      assert(M->Tdt[d]>0);
      M->TDt--;
      heap_pop(&down);
      heap_remove(&up,d);
    } else if ( downv<upv && upv>1.0 ){
      //  increment this
      d = heap_front(&up);
      M->TdT[d]++;
      M->Tdt[d]++;
      assert(M->Tdt[d]<=M->Ndt[d]);
      M->TDt++;
      heap_pop(&up);
      heap_remove(&down,d);
    } else {
      //  none are better so quit
      break;
    }
    if ( M->Tdt[d]<M->Ndt[d] )
      score_up[d] = (ddP.bpar + ddP.apar*M->TdT[d]) 
	/ S_V(ddC.SX,M->Ndt[d],M->Tdt[d]+1);
    else 
      score_up[d] = 0;
    if ( M->Tdt[d]>1 )
      score_down[d] = S_V(ddC.SX,M->Ndt[d],M->Tdt[d])
	/(ddP.bpar + ddP.apar*(M->TdT[d]-1));
    else
      score_down[d] = 0;
    assert(M->Tdt[d]>1||score_down[d]==0);
    assert(M->Tdt[d]<M->Ndt[d] ||score_up[d]==0);
    assert(M->Tdt[d]<=M->Ndt[d]);
    /*
     *  now adjust the two heaps for new vals for [d]
     */
    heap_push(&down,d);
    heap_push(&up,d);
  }  
  free(score_up);
  free(score_down);
  heap_free(&up);
  heap_free(&down);
}
#endif

/*
 *  compute likelihood ratio difference based on *M
 */
static double merge_like(int k1, int k2, merge_t *M) {
  int d;
  double la = 0;
  double lb = log(ddP.bpar);
  int TD_diff = 0;
  double likelihood = 0;
  if ( ddP.apar>0 ) la = log(ddP.apar);
  for (d=0; d<ddN.DT; d++) {
    int Td_diff;  /*  total change in T for doc */
    if ( M->Ndt[d]>1 ) {
      likelihood -= S_S(ddC.SX,ddS.Ndt[d][k2],ddS.Tdt[d][k2]);
      likelihood -= S_S(ddC.SX,ddS.Ndt[d][k1],ddS.Tdt[d][k1]);
      likelihood += S_S(ddC.SX,M->Ndt[d],M->Tdt[d]);
      assert(M->Tdt[d]>=1);
      assert(M->Tdt[d]<=M->Ndt[d]);
      assert(ddS.Ndt[d][k2]==0 || ddS.Tdt[d][k2]>0);
      assert(ddS.Ndt[d][k1]==0 || ddS.Tdt[d][k1]>0);
    }
    yap_infinite(likelihood);
    TD_diff += Td_diff = (M->Tdt[d]-ddS.Tdt[d][k2]-ddS.Tdt[d][k1]);
    if ( Td_diff==0 )
      continue;
    if ( ddP.apar==0 ) {
      likelihood += Td_diff*lb;
    } else {
      likelihood += Td_diff*la;
      if ( Td_diff<0 ) 
	likelihood -= 
	  gammadiff(-Td_diff,M->TdT[d]+ddP.bpar/ddP.apar, 0.0);
      else
	likelihood += 
	  gammadiff(Td_diff, M->TdT[d]-Td_diff+ddP.bpar/ddP.apar, 0.0);
    }
    yap_infinite(likelihood);
  }      
  if ( ddP.PYalpha==H_PDP ) {
    likelihood += (M->TDt-ddS.TDt[k1])*log(ddP.alphapr[k1])
      - ddS.TDt[k2]*log(ddP.alphapr[k2]);
  } else if ( ddP.PYalpha==H_HDP ) {
    likelihood += lgamma(M->TDTm+M->TDt-TD_diff+ddP.b0) 
      - lgamma(M->TDTm+M->TDt+ddP.b0);
    likelihood -= gammadiff(ddS.TDt[k1], ddP.b0*ddP.alphapr[k1], 0.0);
    likelihood -= gammadiff(ddS.TDt[k2], ddP.b0*ddP.alphapr[k2], 0.0);
    likelihood += gammadiff(M->TDt, ddP.b0*ddP.alphapr[k1], 0.0);
  } else {
    double lga0 = lgamma(1-ddP.a0);
    likelihood += lgamma(M->TDTm+M->TDt-TD_diff+ddP.b0) 
      - lgamma(M->TDTm+M->TDt+ddP.b0);
    /*   because k2 gone to zero, so one less topic */
    likelihood -= log(ddP.b0+ddP.a0*(ddS.TDTnz-1));
    if ( ddS.TDt[k2]>1 )
      likelihood -= lgamma(ddS.TDt[k2]-ddP.a0) - lga0;
    if ( ddS.TDt[k1]>1 )
      likelihood -= lgamma(ddS.TDt[k1]-ddP.a0) - lga0;
    likelihood += lgamma(M->TDt-ddP.a0) - lga0;
  }
  yap_infinite(likelihood);
  return likelihood;
}

static double likemerge_alpha(int k1, int k2) {
  /*
   *   storing/saving the new/proposed version
   */
  merge_t M;

  double likelihood;

  if ( k1<=0 || k2<=0 || ddS.TDt[k1]==0 || ddS.TDt[k2]==0 ) 
    return 0.0;
  merge_init(k1, k2, &M);
#ifndef NOOPT_MERGE
  /*  optimise the table counts Tdt[.][k1] */
  merge_opt_Tdt(k1, k2, &M);
#endif
  likelihood = merge_like(k1, k2, &M);
  merge_free(&M);
  return likelihood;
}

typedef struct bestmerge_s {
  int k2;
  double ml;
} bestmerge_t;

static int next_best(bestmerge_t *B) {
  int k, bk=-1;
  double v = 0;
  for (k=0; k<ddN.T; k++) {
    if ( B[k].ml>v ) {
      bk = k;
      v = B[k].ml;
    }
  }
  if ( bk>=0 ) 
    return bk;
  return -1;
}

void like_merge(float minprop, double scale, int best) {
  int k1, k2;
  double likediff;
  float **cmtx;
  int title = 0;
  int mincount = minprop * ddN.NT;
  bestmerge_t B[ddN.T];

  if ( mincount<5 )
    mincount = 5;

  assert(ddP.phi==NULL);
  assert(ddP.theta==NULL);

  for (k1=0; k1<ddN.T; k1++) 
    B[k1].ml = 0;

  cmtx = hca_topmtx();
  if ( !cmtx )
    yap_quit("Out of memory in like_merge()\n");
 
  if ( ddP.PYbeta!=H_None ) 
    yap_quit("Non-parametric beta unimplemented with merge\n");

  for (k1=1; k1<ddN.T; k1++) {
    if ( ddS.NWt[k1]<=mincount ) 
      continue;
    for (k2=0; k2<k1; k2++) {
      // ????? order k1, k2?
      if ( ddS.NWt[k2]<=mincount )
	continue;
      if ( ddP.PYalpha==H_None )
	likediff = likemerge_DIRalpha(k1,k2);
      else
	likediff = likemerge_alpha(k1, k2);
      likediff += likemerge_DIRbeta(k1,k2);
      if ( likediff>0 ) {
	if ( title==0 ) {
	  yap_message("\nPre merge log_2(perp)=%.4lf",  
		      scale * likelihood() );
	}
	if ( verbose>1 ) {
	  if ( title==0 ) 
	    yap_message(", merge report:\n");
	  yap_message("\n   %d+%d cor=%0.6f like+=%0.6g", k1, k2, cmtx[k1][k2], 
		      scale * likediff);
	}     
	title = 1;
	if ( likediff>B[k1].ml ) {
	  B[k1].ml = likediff;
	  B[k1].k2 = k2; 
	}
	if ( likediff>B[k2].ml ) {
	  B[k2].ml = likediff;
	  B[k2].k2 = k1; 
	}
      } else if ( verbose>2 ) {
	  yap_message("\n   %d+%d cor=%0.6f like+=%0.6g", k1, k2, cmtx[k1][k2], 
		      scale * likediff);
      }
 
    }
  }
  while ( best>0 && (k1=next_best(&B[0]))>=0 ) {
    /*
     *   have a good merge at position k1;
     */
    merge_t M;
    yap_message("\n  best merge is %d+%d giving diff=%lf\n", k1, B[k1].k2,
		scale* B[k1].ml);
    merge_init(k1, B[k1].k2, &M);
    hca_merge_stats(k1, B[k1].k2, M.Tdt, NULL);
    // hca_correct_tdt(0);
    merge_free(&M);
    /*  block them from getting picked again */
    B[k1].ml = 0;
    B[B[k1].k2].ml = 0;
  } 
  free(cmtx[0]); free(cmtx);
}
