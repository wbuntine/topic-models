/*
 * Sampling utility for bw
 * Copyright (C) 2011 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *
 *   NOTE:    the auxillary variable method introducing
 *            "q" with a beta distribution means that all b>0,
 *            although a PDP/PDD allows b>-a
 *            i.e.,   our method restricts "b" a little
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

// #define BS_DEBUG

#if 0
/*
 *     these are set in the sample_bw() routine
 *     for use in bterm() and bterm_db()
 */
static  double logsum;

/*
 *   this routine is never called if ddP.aw0==0 && ddP.gamma0==0
 */
static double bw0terms_PDD_old(double bw, void *mydata) {
  double val = -logsum * bw + (PYP_CONC_PSHAPE-1)*log(bw);
  if ( ddP.aw0>0 )
    val += lgamma(bw/ddP.aw0+ddS.TWTnz) - lgamma(bw/ddP.aw0);
  myarms_evals++;
#ifdef BS_DEBUG
  yap_message("Eval bterms(%lf) = %lf\n", bw, val);
#endif
  return val;
}
#endif

static double bw0terms_PDD(double bw, void *mydata) {
  double val = pctl_gammaprior(bw);
  if ( ddP.aw0==0 )
    val += ddS.TWTnz*log(bw);
  else
    val += lgamma(bw/ddP.aw0+ddS.TWTnz) - lgamma(bw/ddP.aw0);
  val -= lgamma(bw+ddS.TWT) - lgamma(bw);
  myarms_evals++;
#ifdef BS_DEBUG
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
#ifdef BS_DEBUG
  yap_message("Eval bterms(%lf) = %lf\n", bw, val);
#endif
  return val;
}

/*
 *    assume a gamma prior of Gamma(1.1,T)
 */
static void sample_bw0_PDD(double *bw) {
  myarmsMH(PYP_CONC_MIN, PYP_CONC_MAX, &bw0terms_PDD, NULL, bw, "bw0", 1);
}

static void sample_bw0_DP(double *bw) {
  /*
   *   compute it in first pass,
   *   then use it inside bterms() and bterms_db()
   */
  assert(ddP.PYbeta==H_HDP);
  myarms(PYP_CONC_MIN, PYP_CONC_MAX, &bw0terms_DP, NULL, bw, "bw0");
}

void sample_bw0(double *bw) {
  double startlike;
  if ( verbose>1 ) {
    startlike = likelihood();
    yap_message("sample_bw0 (pre): bw0=%lf, lp=%lf\n",
                *bw, startlike);
  }
  if ( ddP.PYbeta==H_HDP )
    sample_bw0_DP(&ddP.bw0);
  else 
    sample_bw0_PDD(&ddP.bw0);
#if 0
  if ( *bw<PYP_CONC_MIN )
    *bw = PYP_CONC_MIN;
  if ( *bw>PYP_CONC_MAX )
    *bw = PYP_CONC_MAX;
#endif
  cache_update("bw0");
  if ( verbose>1 ) {
    double endlike = likelihood();
    yap_message("sample_bw0 (post): bw0=%lf, lp=%lf\n",
                *bw, endlike);
    if ( endlike < startlike-50 )
      yap_quit("Sampler failed due to huge decrease!\n");
  }
}
