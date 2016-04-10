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
// #define NOOPT_MERGE

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

//   matches likelihood_DIRalpha()
static double likemerge_DIRalpha(int k1, int k2) {
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

//   likelihood_DIRbeta() 
static double likemerge_DIRbeta(int k1, int k2) {
  int j;
  double val = 0;
  for (j=0; j<ddN.W; j++) {
    if ( ddS.Nwt[j][k2]>0 ) {
      assert(ddP.betapr[j]>0);
      val += (gammadiff((int)ddS.Nwt[j][k1]+ddS.Nwt[j][k2], ddP.betapr[j], 0.0)
	      - gammadiff((int)ddS.Nwt[j][k1], ddP.betapr[j], 0.0)
	      - gammadiff((int)ddS.Nwt[j][k2], ddP.betapr[j], 0.0));
    } 
  }     
  val -= (gammadiff((int)ddS.NWt[k1]+ddS.NWt[k2], ddP.betatot, 0.0)
	  - gammadiff((int)ddS.NWt[k1], ddP.betatot, 0.0)
	  - gammadiff((int)ddS.NWt[k2], ddP.betatot, 0.0));
  yap_infinite(val);
  return val;
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
typedef struct merge_alpha_s {
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
} merge_alpha_t;

typedef struct merge_beta_s {
  /*
   *   storing/saving the new/proposed version
   */
  uint32_t NWt, TWt, TWTm;
  /*
   *  cache for doc side change
   */
  uint32_t *TwT;
  uint32_t *Nwt;
  uint16_t *Twt;
} merge_beta_t;


static void merge_free_Tdt(merge_alpha_t *M) {
  free(M->Ndt);
  free(M->TdT);
  free(M->Tdt);
}
static void merge_free_Twt(merge_beta_t *M) {
  free(M->Nwt);
  free(M->TwT);
  free(M->Twt);
}

static void merge_init_Tdt(int k1, int k2, merge_alpha_t *M) {
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
  M->TDt = ddS.TDt[k1] + ddS.TDt[k2];
  for (d=0; d<ddN.DT; d++) {
    M->Tdt[d] = ddS.Tdt[d][k1] + ddS.Tdt[d][k2];
    M->Ndt[d] = ddS.Ndt[d][k1] + ddS.Ndt[d][k2];
    M->TdT[d] = 0;
    for (k=0; k<ddN.T; k++)
      M->TdT[d] += ddS.Tdt[d][k];
    assert(M->Tdt[d]<=M->Ndt[d]);
  }
  M->TDTm = ddS.TDT - M->TDt;
#ifndef NDEBUG
  for (d=0; d<ddN.DT; d++) {
    assert(M->Tdt[d]<=M->Ndt[d]);
  }
#endif
}

static void merge_init_Twt(int k1, int k2, merge_beta_t *M) {
  int w;
  /*
   *   build local store
   */
  M->Twt = u16vec(ddN.W); 
  M->TwT = u32vec(ddN.W); 
  M->Nwt = u32vec(ddN.W);
  if ( !M->Nwt || !M->TwT || !M->Twt  )
    yap_quit("Out of memory in likemerge()\n");

  /*
   *  initialise all data entries to simple sum
   */
  M->NWt = ddS.NWt[k1] + ddS.NWt[k2];
  M->TWt = ddS.TWt[k1] + ddS.TWt[k2];
  for (w=0; w<ddN.W; w++) {
    M->Twt[w] = ddS.Twt[w][k1] + ddS.Twt[w][k2];
    M->Nwt[w] = ddS.Nwt[w][k1] + ddS.Nwt[w][k2];
    M->TwT[w] = ddS.TwT[w];
  }
  M->TWTm = ddS.TWT - M->TWt;
#ifndef NDEBUG
  for (w=0; w<ddN.W; w++) {
    assert(M->Twt[w]<=M->Nwt[w]);
  }
#endif
}

#ifdef UNMADE
#define K_BOUND 10
static int getTdT(int d) {
  int k1;
  int TdT=0;
  uint16_t *Tdt = ddS.Tdt[d];
  for (k=0; k<ddN.T; k++)
    TdT += Tdt[k];
  return TdT;
}
static double merge_sumapprox_Tdt(int k) {
  double x[ddN.DT];
  int d;
  for (d=0; d<ddN.DT; d++) {
    double xp = 0, xm = 0;
    int nn, tt, TdT;
    nn = ddS.Ndt[d][k];
    if ( nn<=1 ) {
      x[d] = 0;
      continue;
    }
    tt = ddS.Tdt[d][k];
    TdT = getTdT(d);
    if ( tt<nn )
      xp = (ddP.bpar + ddP.apar*TdT) * S_V(ddC.SX,nn,tt+1)
	* alphabasetopicprob(k);
    if ( tt>1) {
      ddS.TDt[k]--; ddS.TDT--;
      xm = 1.0/(ddP.bpar + ddP.apar*(TdT-1)) / S_V(ddC.SX,nn,tt)
	/ alphabasetopicprob(k);
      ddS.TDt[k]++; ddS.TDT++;
    }
    if ( xm>xp )
      x[d] = xm;
    else 
      x[d] = xp;
  }
  // ?????????
}
#endif

#ifndef NOOPT_MERGE
static void merge_opt_Tdt(int k1, int k2, merge_alpha_t *M) {
int d;
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
	* S_V(ddC.SX,M->Ndt[d],M->Tdt[d]+1);
    else 
      score_up[d] = 0;
    if ( M->Tdt[d]>1 )
      score_down[d] = 1.0 / S_V(ddC.SX,M->Ndt[d],M->Tdt[d])
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
	* S_V(ddC.SX,M->Ndt[d],M->Tdt[d]+1);
    else 
      score_up[d] = 0;
    if ( M->Tdt[d]>1 )
      score_down[d] = 1.0 / S_V(ddC.SX,M->Ndt[d],M->Tdt[d])
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
static double merge_like_Twt(int k1, int k2, merge_beta_t *M) {
  int i, w;
  double likelihood = 0;
#ifndef BWPAR0
  double lbw = log(ddP.bwpar);
#endif
  double law = log(ddP.awpar);
  double TW_diff = 0;
#ifdef BWPAR0
	yap_quit("BWPAR0 unimpleented in merge\n");
#endif
  for (i=0; i<ddN.W; i++) {
    likelihood -= S_S(ddC.SY,ddS.Nwt[i][k1],ddS.Twt[i][k2]);
    likelihood -= S_S(ddC.SY,ddS.Nwt[i][k1],ddS.Twt[i][k2]);
    likelihood += S_S(ddC.SY,M->Nwt[i],M->Twt[i]);
  }
  if ( ddP.awpar==0 ) {
#ifdef BWPAR0
    likelihood += M->TWt-ddS.TWt[k1]*log(ddP_bwpar(k1))
	-ddS.TWt[k2]*log(ddP_bwpar(k2));
#else
    likelihood += (M->TWt-ddS.TWt[k1]-ddS.TWt[k2])*lbw;
#endif
  } else {
    likelihood += (M->TWt-ddS.TWt[k1]-ddS.TWt[k2])*law 
      + gammadiff((int)M->TWt, ddP.bwpar/ddP.awpar, 0.0)
      - gammadiff((int)ddS.TWt[k1], ddP_bwpar(k1)/ddP.awpar, 0.0)
      - gammadiff((int)ddS.TWt[k2], ddP_bwpar(k2)/ddP.awpar, 0.0);
  }
  likelihood += gammadiff((int)ddS.NWt[k1], ddP_bwpar(k1), 0.0);
  likelihood += gammadiff((int)ddS.NWt[k2], ddP_bwpar(k2), 0.0);
  likelihood -= gammadiff((int)M->NWt, ddP.bwpar, 0.0);
  yap_infinite(likelihood);
  if ( ddP.PYbeta==H_PDP ) {
    for (w=0; w<ddN.W; w++) {
      if ( ddS.TwT[w]>0 ) {
	// ???????????????
        likelihood += ddS.TwT[w]*log(ddP.betapr[w]);
      }
    }      
  } else if ( ddP.PYbeta==H_HDP ) {
    yap_quit("merge with PYbeta unimplemented\n");
    likelihood += lgamma(M->TWTm+M->TWt-TW_diff+ddP.bw0) 
      - lgamma(M->TWTm+M->TWt+ddP.bw0);
    for (w=0; w<ddN.W; w++) {
      // ???????????
      likelihood -= gammadiff(ddS.TWt[k1], ddP.bw0*ddP.betapr[k1], 0.0);
      likelihood -= gammadiff(ddS.TWt[k2], ddP.bw0*ddP.betapr[k2], 0.0);
      likelihood += gammadiff(M->TWt, ddP.bw0*ddP.betapr[k1], 0.0);
    }
  } else {
    double lgaw0 = lgamma(1-ddP.aw0);
    likelihood += lgamma(M->TWTm+M->TWt-TW_diff+ddP.bw0) 
      - lgamma(M->TWTm+M->TWt+ddP.bw0);
    /*   because k2 gone to zero, so one less topic */
    likelihood -= log(ddP.bw0+ddP.aw0*(ddS.TWTnz-1));
    if ( ddS.TWt[k2]>1 )
      likelihood -= lgamma(ddS.TWt[k2]-ddP.aw0) - lgaw0;
    if ( ddS.TWt[k1]>1 )
      likelihood -= lgamma(ddS.TWt[k1]-ddP.aw0) - lgaw0;
    likelihood += lgamma(M->TWt-ddP.aw0) - lgaw0;
  }
  yap_infinite(likelihood);
  return likelihood;
}

#if 0
/*
 *  compute likelihood ratio difference based on *M
 */
static double merge_like_Tdt_sum(int k1, int k2, merge_alpha_t *M) {
  double *val;
  int d;
  val = dvec(ddN.DT);

  for (d=0; d<ddN.DT; d++) {

}
#endif

static double merge_like_Tdt(int k1, int k2, merge_alpha_t *M) {
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
  merge_alpha_t M;
  double likelihood;

  if ( k1<=0 || k2<=0 || ddS.TDt[k1]==0 || ddS.TDt[k2]==0 ) 
    return 0.0;
  merge_init_Tdt(k1, k2, &M);
#ifndef NOOPT_MERGE
  /*  optimise the table counts Tdt[.][k1] */
  merge_opt_Tdt(k1, k2, &M);
#endif
  likelihood = merge_like_Tdt(k1, k2, &M);
  merge_free_Tdt(&M);
  return likelihood;
}

static double likemerge_beta(int k1, int k2) {
  /*
   *   storing/saving the new/proposed version
   */
  merge_beta_t M;

  double likelihood;

  if ( k1<=0 || k2<=0 || ddS.TWt[k1]==0 || ddS.TWt[k2]==0 ) 
    return 0.0;
  merge_init_Twt(k1, k2, &M);
#ifndef NOOPT_MERGE
  /*  optimise the table counts Tdt[.][k1] */
  //  ???????merge_opt_Twt(k1, k2, &M);
#endif
  likelihood = merge_like_Twt(k1, k2, &M);
  merge_free_Twt(&M);
  return likelihood;
}

typedef struct bestmerge_s {
  int k2;
  double ml;
} bestmerge_t;

static int next_best(bestmerge_t *B) {
  int k, bk=-1;
  double v = 0;
  yap_message("merge buffer: ");
  for (k=0; k<ddN.T; k++) {
    if ( B[k].ml>0 ) {
      yap_message("%d+%d ", k, B[k].k2);
    }
  }
  yap_message("\n");
  for (k=0; k<ddN.T; k++) {
    if ( B[k].ml>v ) {
      bk = k;
      v = B[k].ml;
    }
  }
  return bk;
}

void like_merge(float minprop, double scale, int best) {
  int k1, k2;
  double realdiff = 0;
  double likediff;
  int got=0;
  /*  only use this for reporting ; should disable in production */
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
      if ( ddS.NWt[k2]<=mincount )
	continue;
      /*  now have a pair to check (k1,k2) with OK counts */
      if ( ddP.PYalpha==H_None )
	likediff = likemerge_DIRalpha(k1,k2);
      else
	likediff = likemerge_alpha(k1, k2);
      if ( ddP.PYbeta==H_None )
        likediff += likemerge_DIRbeta(k1,k2);
      else
	likediff += likemerge_beta(k1, k2);
      if ( likediff>0 ) {
	got++;
	if ( title==0 && verbose ) {
	  double like = scale * likelihood();
	  yap_message("\nPre merge log_2(perp)=%.4lf", like);
	  realdiff = like;
	}
	if ( verbose>1 ) {
	  if ( title==0 ) 
	    yap_message(", merge report:\n");
	  yap_message("   %d+%d cor=%0.6f like+=%0.6g", k1, k2, cmtx[k1][k2], 
		      scale * likediff);
	}     
	title = 1;
	if ( likediff>B[k1].ml ) {
	  B[k1].ml = likediff;
	  B[k1].k2 = k2; 
	  if ( verbose>1 ) yap_message("  stored");
	}
	if ( likediff>B[k2].ml ) {
	  B[k2].ml = likediff;
	  B[k2].k2 = k1; 
	  if ( verbose>1 ) yap_message("  stored");
	}
	if ( verbose>1 ) yap_message("\n");
      } else if ( verbose>2 ) {
        yap_message("   %d+%d cor=%0.6f like+=%0.6g\n", k1, k2, cmtx[k1][k2], 
                    scale * likediff);
      }
    }
  }
  while ( got && best-->0 && (k1=next_best(&B[0]))>=0 ) {
    /*
     *   have a good merge at position k1;
     */
    merge_alpha_t Ma;
    merge_beta_t Mb;
    k2 = B[k1].k2;
    yap_message("  best merge is %d+%d giving diff=%lf\n", k1, k2,
		scale* B[k1].ml);
    Ma.Tdt = NULL;
    Mb.Twt = NULL;
    if ( ddP.PYalpha ) 
      merge_init_Tdt(k1, k2, &Ma);
    if ( ddP.PYbeta ) 
      merge_init_Twt(k1, k2, &Mb);	       
    //  WRAY:  need to checkk what this does, it it why change?
    hca_merge_stats(k1, k2, Ma.Tdt, Mb.Twt);
    // correct_tdt(0);
    if ( ddP.PYalpha ) 
      merge_free_Tdt(&Ma);
    if ( ddP.PYbeta ) 
      merge_free_Twt(&Mb);
    /*  block them from getting picked again */
    B[k1].ml = 0;
    B[k2].ml = 0;
    {
      int k;
      for (k=0; k<ddN.T; k++) {
	if ( B[k].k2==k1 || B[k].k2==k2 )
	  B[k].ml = 0;
      } 
    }
  }	
  if ( got && verbose ) {
    double like = scale * likelihood();
    realdiff -= like;
    yap_message("\nPost merge log_2(perp)=%.4lf (%.6lf)", like, realdiff);
  }  
  if ( got==0 && verbose ) {
    yap_message("Merge found no candidates\n");
  }
  free(cmtx[0]); free(cmtx);
}
