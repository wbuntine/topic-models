/*
 * Sampling utility for b
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
#include "sample.h"
#include "hca.h"
#include "stats.h"


// #define B_DEBUG


static double b0terms_PDD(double b, void *mydata) {
  double val = pctl_gammaprior(b);
  if ( ddP.a0==0 )
    val += ddS.TDTnz*log(b);
  else
    val += lgamma(b/ddP.a0+ddS.TDTnz) - lgamma(b/ddP.a0);
  val -= lgamma(b+ddS.TDT) - lgamma(b);
  myarms_evals++;
#ifdef B_DEBUG
  yap_message("Eval b0terms_PDD(%lf) = %lf\n", b, val);
#endif
  return val;
}

static double b0terms_DP(double b, void *mydata) {
  double val = - log(b);
  int t;
  for (t=0; t<ddN.T; t++)
    if ( ddS.TDt[t]>0 ) {
      val += gammadiff((int)ddS.TDt[t], b/ddN.T, 0.0);
    }
  val -= lgamma(b+ddS.TDT) - lgamma(b);
  myarms_evals++;
#ifdef B_DEBUG
  yap_message("Eval b0terms_DP(%lf) = %lf\n", b, val);
#endif
  return val;
}

/*
 *
 */
static double bterms(double b, void *mydata) {
  int i;
  uint16_t *localTd = (uint16_t *)mydata;
  double val = pctl_gammaprior(b);
  double lgb = lgamma(b);
  double lgba = 0;
  if ( ddP.apar>0 )
    lgba = lgamma(b/ddP.apar);
  for (i=0; i<ddN.DT; i++) {
    if ( ddP.apar>0 )
      val += gammadiff(localTd[i], b/ddP.apar, lgba);
    else 
      val += localTd[i] * log(b);
    val -= gammadiff((int)ddS.NdT[i], b, lgb) ;
  }
  myarms_evals++;
#ifdef B_DEBUG
  yap_message("Eval bterms(%lf) = %lf", b, val);
  ddP.bpar = b;
  cache_update("b");
  yap_message(", lp=%lf\n", likelihood());
#endif
  return val;
}

/*
 *    just call the likelihood function
 */
struct bdkterms_s {
  int t;
  uint16_t **docstats;
};
static double bdkterms(double b, void *mydata) {
  struct bdkterms_s *ps = (struct bdkterms_s*)mydata;
  int t = ps->t;
  double val;
  ddP.bdk[t] = b;
  val = dmi_likelihood_bterms(&ddM, t, ps->docstats, pctl_gammaprior,
                              ddP.ad, ddP.bdk);
  myarms_evals++;
#ifdef B_DEBUG
  yap_message("Eval (from likelihood) bdk[%d]terms(%lf) = %lf", t,b,val);
#endif
  return val;
}

/*
 */
static double betaterms(double mytbeta, void *mydata) {
  int j,t;
  double val = 0;
  double old_beta = *(double*)mydata;
#ifdef B_DEBUG
  double like;
#endif
  for (t=0; t<ddN.T; t++) {
    for (j=0; j<ddN.W; j++) {
      if ( ddS.Nwt[j][t]>0 ) {
        val += gammadiff((int)ddS.Nwt[j][t], mytbeta*ddP.betapr[j]/old_beta, 0);
      }
    }
    val -= gammadiff((int)ddS.NWt[t], mytbeta, 0);
  }
  myarms_evals++;
  myarms_last = mytbeta;
#ifdef B_DEBUG
  yap_message("Eval betaterms(%lf) = %lf", mytbeta, val);
  ddP.beta = mytbeta;
  cache_update("beta");
  old_beta = mytbeta;
  like = likelihood();
  if ( last_val != 0 ) {
    yap_message(", lp=%lf diffs=%lf vs %lf\n", like,
                val-last_val, like-last_like);
  }
  last_like = like;
  last_val = val;
#endif
  return val;
}

static double bw0terms_PDD(double bw, void *mydata) {
  double val = pctl_gammaprior(bw);
  if ( ddP.aw0==0 )
    val += ddS.TWTnz*log(bw);
  else
    val += lgamma(bw/ddP.aw0+ddS.TWTnz) - lgamma(bw/ddP.aw0);
  val -= lgamma(bw+ddS.TWT) - lgamma(bw);
  myarms_evals++;
#ifdef B_DEBUG
  yap_message("Eval bterms(%lf) = %lf\n", bw, val);
#endif
  return val;
}

static double bw0terms_DP(double bw, void *mydata) {
  double val = - log(bw);
  int j;
  for (j=0; j<ddN.W; j++)
    if ( ddS.TwT[j]>0 ) {
      val += gammadiff((int)ddS.TwT[j], bw*ddP.betapr[j], 0.0);
    }
  val -= lgamma(bw+ddS.TWT) - lgamma(bw);
  myarms_evals++;
#ifdef B_DEBUG
  yap_message("Eval bterms(%lf) = %lf\n", bw, val);
#endif
  return val;
}

static double bwterms(double bw, void *mydata) {
  int t;
  double val = pctl_gammaprior(bw);
#ifdef SBW_USECACHE
  struct gcache_s lgba_t;
  struct gcache_s lgb_t;
#else
  double lgb = lgamma(bw);
  double lgba = 0;
#endif
#ifdef SBW_USECACHE
  if ( ddP.awpar>0 )
    gcache_init(&lgba_t, bw/ddP.awpar);
  gcache_init(&lgb_t, bw);
#else
  if ( ddP.awpar>0 )
    lgba = lgamma(bw/ddP.awpar);
#endif
  for (t=0; t<ddN.T; t++) {
#ifdef SBW_USECACHE
    val += gcache_value(&lgba_t, (int)ddS.TWt[t])
      - gcache_value(&lgb_t, (int)ddS.NWt[t]);
#else
    if ( ddP.awpar>0 )
      val += gammadiff(ddS.TWt[t], bw/ddP.awpar, lgba);
    else
      val += ddS.TWt[t] * log(bw);
    val -= gammadiff((int)ddS.NWt[t], bw, lgb) ;
#endif
  }
  myarms_evals++;
#ifdef B_DEBUG
  yap_message("Eval bwterms(%lf) = %lf", bw, val);
  ddP.bwpar = bw;
  cache_update("bw");
  yap_message(", lp=%lf\n", likelihood());
#endif
  return val;
}

/************************************************************
 *
 *   main routines
 *
 ************************************************************/

void sample_b0(double *b) {
  if ( ddP.PYbeta==H_HDP ) {
    assert(ddP.PYbeta==H_HDP);
    myarmsMH(PYP_CONC_MIN, PYP_CONC_MAX, &b0terms_DP, NULL, b, "b0", 1);
  } else {
    /*
     *  prior is pctl_gammaprior
     */
    myarmsMH(PYP_CONC_MIN, PYP_CONC_MAX, &b0terms_PDD, NULL, b, "b0", 1);
  }
  cache_update("b0");
}

/*
 *    this is the sampler given in Lan Du's papers
 */
void sample_b(double *b) {
  int i, t;  
  uint16_t *localTd; 
  /*
   *   compute t totals for docs since not stored
   */
  localTd = malloc(sizeof(*localTd)*ddN.DT);
  for (i=0; i<ddN.DT; i++) {
    uint16_t Td_ = 0;
    for (t=0; t<ddN.T; t++)
      Td_ += ddS.Tdt[i][t];
    localTd[i] = Td_;
  }
  
  myarmsMH(PYP_CONC_MIN, PYP_CONC_MAX, &bterms, localTd, b, "b", 1);
  cache_update("b");
  free(localTd);
}

void sample_bdk(double *b, int k) {
  struct bdkterms_s ps;
  ps.t = k;
  assert(b);
  ps.docstats = ddP.docstats;
  myarmsMH(PYP_CONC_MIN, PYP_CONC_MAX, &bdkterms, &ps, &b[k], "bdk", 1);
}

/*
 *  this allows Dirichlet prior to be non-uniform,
 *  so optimisation done on total weight;
 */
void sample_beta(double *mytbeta) {
  double bmax =  DIR_MAX*ddN.W;
  double old_beta = ddP.beta;
#ifdef B_DEBUG
  last_val = 0;
  last_like = 0;
#endif
  // if ( bmax > DIR_TOTAL_MAX )  bmax = DIR_TOTAL_MAX;
  /*
   *  beta is set using ddP.beta, so need to factor it out
   */
  myarmsMH(DIR_MIN*ddN.W, bmax, &betaterms, &old_beta, mytbeta, "beta",1);
  cache_update("beta");
}


void sample_bw0(double *bw) {
  if ( ddP.PYbeta==H_HDP ) {
    assert(ddP.PYbeta==H_HDP);
    myarms(PYP_CONC_MIN, PYP_CONC_MAX, &bw0terms_DP, NULL, bw, "bw0");
  } else {
    /*
     *    assume a gamma prior pctl_gammaprior()
     */
    myarmsMH(PYP_CONC_MIN, PYP_CONC_MAX, &bw0terms_PDD, NULL, bw, "bw0", 1);
  }
  cache_update("bw0");
}

void sample_bw(double *bw) {
  myarmsMH(PYP_CONC_MIN, PYP_CONC_MAX, &bwterms, NULL, bw, "bw", 1);
  cache_update("bw");
}
