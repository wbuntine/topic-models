/*
 * Sampling utility for a
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

// #define A_DEBUG
// #define A1_DEBUG

#ifdef A1_DEBUG
  static double last_val = 0;
  static double last_like = 0;
#endif
/*
 *  globals from hca.h and stats.h
 */
extern int verbose;



/*
 */
static double a0terms(double mya0, void *mydata) {
  int i;
  double l1a0 = log(1-mya0);
  double l2a0 = log((1-mya0)*(2-mya0));
  double lga0 = lgamma(1-mya0);
  double val = 0;
#ifdef A_DEBUG
  double like;
#endif
  val += ddS.TDTnz*log(mya0) + lgamma(ddP.b0/mya0+ddS.TDTnz)
    - lgamma(ddP.b0/mya0);
  for (i=0; i<ddN.T; i++)
    /*  note the root node is a PDD so all t's = 1 */
    if ( ddS.TDt[i]>1 ) {
      if ( ddS.TDt[i]==2 )
        val += l1a0;
      else if ( ddS.TDt[i]==3 )
        val += l2a0;
      else
        val += lgamma(ddS.TDt[i]-mya0) - lga0;
    }
#ifdef A_DEBUG
  yap_message("Eval a0terms(%lf) = %lf", mya0, val);
  ddP.a0 = mya0;
  cache_update("a0");
  like = likelihood();
  if ( last_val != 0 ) {
    yap_message(", lp=%lf diffs=%lf vs %lf\n", like,
                val-last_val, like-last_like);
  }
  last_like = like;
  last_val = val;
#endif
  myarms_evals++;
  return val;
}

/*
 */
static double aterms(double mya, void *mydata) {
  int i, t;
  double val = 0;
  double la = log(mya);
#ifdef A_DEBUG
  float save_a = ddC.SX->a;
  double like;
#endif
  S_remake(ddC.SX, mya);
  for (i=0; i<ddN.DT; i++) {
    uint32_t Td_ = 0;
    for (t=0; t<ddN.T; t++) {
      Td_ += ddS.Tdt[i][t];
      if ( ddS.Ndt[i][t]>1 ) {
	val += S_S(ddC.SX,ddS.Ndt[i][t],ddS.Tdt[i][t]);
      }
    }
    val += Td_*la + lgamma(ddP.bpar/mya+Td_) - lgamma(ddP.bpar/mya);
  }  
  myarms_evals++;
#ifdef A_DEBUG
  yap_message("Eval aterms(%lf) = %lf (S had %f)", mya, val, save_a);
  ddP.apar = mya;
  cache_update("a");
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

/*
 */
static double adkterms(double mya, void *mydata) {
  double val = 0;
  uint16_t **docstats = (uint16_t **)mydata;
#ifdef A_DEBUG
  float save_a = ddC.SD->a;
  double like;
#endif
  ddP.ad = mya;
  cache_update("ad");
  val = dmi_likelihood_aterms(&ddM, docstats,
                              pctl_gammaprior, ddP.ad, ddP.bdk, ddC.SD);
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

#define CONJPRIOR
static double alphaterms(double alphatot, void *mydata) {
  int t,s;
  double val = 0;
  double tot;
  double lga = lgamma(alphatot/ddN.T);
  double lgat = lgamma(alphatot); 
#ifdef A_DEBUG
  double like;
#endif
#ifdef CONJPRIOR
  val += ddN.T*(lgamma((alphatot+1.0)/ddN.T) - lga);
  val -= lgamma(alphatot+1.0) - lgamma(alphatot);
#endif
  for (s=0; s<ddN.DT; s++) {
    tot = 0;
    for (t=0; t<ddN.T; t++) {
      tot += alphatot/ddN.T+ddS.Ndt[s][t];
      val += gammadiff(ddS.Ndt[s][t],alphatot/ddN.T,lga);
    }
    val -= lgamma(tot) - lgat;
  }
  myarms_evals++;
  myarms_last = alphatot;
#ifdef A_DEBUG
  yap_message("Eval alphaterms(%lf) = %lf\n", alphatot, val);
  ddP.alphatot = alphatot;
  cache_update("alpha");
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

/*  
 */
static double aw0terms(double myaw0, void *mydata) {
  int i;
  double l1aw0 = log(1-myaw0);
  double l2aw0 = log((1-myaw0)*(2-myaw0));
  double lgaw0 = lgamma(1-myaw0);
  double val = 0;
#ifdef A_DEBUG
  double like;
#endif
  val += ddS.TWTnz*log(myaw0) + lgamma(ddP.bw0/myaw0+ddS.TWTnz)
    - lgamma(ddP.bw0/myaw0);
  for (i=0; i<ddN.W; i++)
    /*  note the root node is a PDD so all t's = 1 */
    if ( ddS.TwT[i]>1 ) {
      if ( ddS.TwT[i]==2 )
        val += l1aw0;
      else if ( ddS.TwT[i]==3 )
        val += l2aw0;
      else
        val += lgamma(ddS.TwT[i]-myaw0) - lgaw0;
    }
#ifdef A_DEBUG
  yap_message("Eval aw0terms(%lf) = %lf", myaw0, val);
  ddP.aw0 = myaw0;
  cache_update("aw0");
  like = likelihood();
  if ( last_val != 0 ) {
    yap_message(", lp=%lf diffs=%lf vs %lf\n", like,
                val-last_val, like-last_like);
  }
  last_like = like;
  last_val = val;
#endif
  myarms_evals++;
  return val;
}

/************************************************************
 *
 *   main routines
 *
 ************************************************************/

void sample_a0(double *mya0) {
#ifdef A_DEBUG
  last_val = 0;
  last_like = 0;
#endif
  myarms(PYP_DISC_MIN, PYP_DISC_MAX, &a0terms, NULL, mya0, "a0");
  cache_update("a0");
}

void sample_a(double *mya) {
#ifdef A_DEBUG
  last_val = 0;
  last_like = 0;
#endif
  if ( verbose>1 )
    yap_message("sample_a (pre):  a=%lf, lp=%lf\n",
		*mya, likelihood());
  myarms(PYP_DISC_MIN, PYP_DISC_MAX, &aterms, NULL, mya, "a");
  cache_update("a");
  if ( verbose>1 )
    yap_message("sample_a (post):  a=%lf, lp=%lf\n",
		*mya, likelihood());
}

void sample_adk(double *mya) {
  uint16_t **docstats;
#ifdef A_DEBUG
  last_val = 0;
  last_like = 0;
#endif
  docstats = dmi_astore(&ddM);
  myarms(PYP_DISC_MIN, PYP_DISC_MAX, &adkterms, docstats, mya, "adk");
  cache_update("ad");
  dmi_freeastore(&ddM, docstats);
}

/*
 *  assumes uniform prior Dirichlet
 */
void sample_alpha(double *alphatot) {
  double dirmax = DIR_TOTAL_MAX;
  if ( dirmax>ddN.T * DIR_MAX )
    dirmax = ddN.T * DIR_MAX;
#ifdef A_DEBUG
  last_val = 0;
  last_like = 0;
#endif
  if ( myarmsMH(DIR_MIN*ddN.T, dirmax,
                &alphaterms, NULL, alphatot, "alphatot",1) ) {
    yap_message("sample_alpha: error in result\n");
  }
  cache_update("alpha");
}

void sample_aw0(double *myaw0) {
#ifdef A_DEBUG
  last_val = 0;
  last_like = 0;
#endif
  /*
   *   compute it in first pass,
   *   then use it inside aw0terms() and aw0terms_da()
   */
  myarms(PYP_DISC_MIN, PYP_DISC_MAX, &aw0terms, NULL, myaw0, "aw0");
  cache_update("aw0");
}


/*
 */
static double awterms(double myaw, void *mydata) {
  int i, t;
  double val = 0;
  double law = log(myaw);
#ifdef A_DEBUG
  float save_a = ddC.SY->a;
  double like;
#endif
  S_remake(ddC.SY, myaw);
  for (t=0; t<ddN.T; t++) {
    uint32_t Tw_ = 0;
    for (i=0; i<ddN.W; i++) {
      Tw_ += ddS.Twt[i][t];
      if ( ddS.Nwt[i][t]>1 ) {
        val += S_S(ddC.SY,ddS.Nwt[i][t],ddS.Twt[i][t]);
      }
    }
    val += Tw_*law + lgamma(ddP_bwpar(t)/myaw+Tw_) - lgamma(ddP_bwpar(t)/myaw);
  }
  myarms_evals++;
#ifdef A_DEBUG
  yap_message("Eval awterms(%lf) = %lf (S had %f)", myaw, val, save_a);
  ddP.awpar = myaw;
  cache_update("aw");
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

void sample_aw(double *myaw) {
#ifdef A_DEBUG
  last_val = 0;
  last_like = 0;
#endif    
  /*
   *   compute it in first pass,
   *   then use it inside awterms() and awterms_da()
   */
  myarms(PYP_DISC_MIN, PYP_DISC_MAX, &awterms, NULL, myaw, "aw");
  cache_update("aw");
}


/*************************************************/

static double ngashterms(double myng, void *mydata) {
  int t;
  double val;
#ifdef A1_DEBUG
  double like;
#endif
  val = PYP_CONC_GAMMA(myng);
  for (t=0; t<ddN.T; t++) {
    val += lgamma(myng+ddS.TDt[t]) - lgamma(myng);
    val -= myng*log(1.0+ddS.NGscalestats[t]*ddP.ngasc);
  }
#ifdef A1_DEBUG
  yap_message("Eval ngashterms(%lf) = %lf", myng, val);
  ddP.ngash = myng;
  like = likelihood();
  if ( last_val != 0 ) {
    yap_message(", lp=%lf diffs=%lf vs %lf\n", like,
                val-last_val, like-last_like);
}
  last_like = like;
  last_val = val;
#endif
  myarms_evals++;
  return val;
}
  
static double ngascterms(double myng, void *mydata) {
  int t;
  double val;
#ifdef A1_DEBUG
  double like;
#endif
  val = PYP_CONC_GAMMA(1.0/myng);
  for (t=0; t<ddN.T; t++) {
    val -= ddP.ngash*log(myng);
    val -= (ddP.ngash+ddS.TDt[t]) *
      log(1.0/myng + ddS.NGscalestats[t]);
  }
#ifdef A1_DEBUG
  yap_message("Eval ngascterms(%lf) = %lf", myng, val);
  ddP.ngasc = myng;
  like = likelihood();
  if ( last_val != 0 ) {
    yap_message(", lp=%lf diffs=%lf vs %lf\n", like,
                val-last_val, like-last_like);
  }
  last_like = like;
  last_val = val;
#endif
  myarms_evals++;
  return val;
}

void sample_ngash(double *myng) {
#ifdef A1_DEBUG
  last_val = 0;
  last_like = 0;
#endif
  NGscalestats(1);
  myarms(PYP_CONC_MIN, PYP_CONC_MAX, &ngashterms, NULL, myng, "ngash");
  ddP.ngash = *myng;
  cache_update("ngash");
}
void sample_ngasc(double *myng) {
#ifdef A1_DEBUG
  last_val = 0;
  last_like = 0;
#endif
  NGscalestats(1);
  myarmsMH(0.01, 100, &ngascterms, NULL, myng, "ngasc", 1);
  ddP.ngasc = *myng;
  cache_update("ngasc");
}


/*****************************************************/


static double UNterms(double myUN, void *mydata) {
  int d = *(int*)mydata;
  int k;
  double like = 0;
  like = (ddS.NdT[d]-1) * log(myUN);
  for (k=0; k<ddN.T; k++)
#ifdef NG_SPARSE
    if (  M_docsparse(d,k) )
#endif
      like -= (ddP.alphapr[k]+ddS.Ndt[d][k])*log(ddS.UN[d]+ddP.NGbeta[k]);
  return like;
}

void sample_UN(int d) {
  int k;
  /*
   *   compute it in first pass,
   *   then use it inside awterms() and awterms_da()
   */
  for (k=0; k<ddN.T; k++)
    ddP.alphapr[k] = (ddP.ngash+ddS.TDt[k])/(1/ddP.ngasc+ddS.NGscalestats[k]);
  myarms(0.00001, ddN.NT,  &UNterms, (void*)&d, &ddS.UN[d], "UN");
}

#define UN_SAMPLE
void opt_UN(int did) {
    double val = 0;
    int t;
    assert(ddS.UN);
    for (t=0; t<ddN.T; t++) {
      double alphaprt;
#ifdef NG_SPARSE
      if (  M_docsparse(did,t) ) 
#endif
	alphaprt = (ddP.ngash+ddS.TDt[t])/(1/ddP.ngasc+ddS.NGscalestats[t]);
	val += ((double)ddS.Ndt[did][t]+alphaprt)
	  / (ddS.UN[did]+ddP.NGbeta[t]);
    }
    assert(val>0);
#ifdef UN_SAMPLE
    //   approximate sampler
    ddS.UN[did] = rng_gamma(rngp,ddS.NdT[did])/val;
#else
    //   fixed point
    ddS.UN[did] = (ddS.NdT[did]-1.0)/val;
#endif
    assert(ddS.UN[did]>0);
}
