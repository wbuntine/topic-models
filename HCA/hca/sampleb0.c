/*
 * Sampling utility for b0
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

// #define BS_DEBUG

/*
 *  globals from hca.h and stats.h
 */

static double b0terms_PDD(double b, void *mydata) {
  double val = pctl_gammaprior(b);
  if ( ddP.a0==0 )
    val += ddS.TDTnz*log(b);
  else
    val += lgamma(b/ddP.a0+ddS.TDTnz) - lgamma(b/ddP.a0);
  val -= lgamma(b+ddS.TDT) - lgamma(b);
  myarms_evals++;
#ifdef BS_DEBUG
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
#ifdef BS_DEBUG
  yap_message("Eval b0terms_DP(%lf) = %lf\n", b, val);
#endif
  return val;
}

/*
 *    assume a gamma prior of Gamma(1.1,T)
 */
static void sample_b0_PDD(double *b) {
  myarmsMH(PYP_CONC_MIN, PYP_CONC_MAX, &b0terms_PDD, NULL, b, "b0", 1);
}

static void sample_b0_DP(double *b) {
  /*
   *   compute it in first pass,
   *   then use it inside bterms() and bterms_db()
   */
  assert(ddP.PYbeta==H_HDP);
  myarmsMH(PYP_CONC_MIN, PYP_CONC_MAX, &b0terms_DP, NULL, b, "b0", 1);
}

void sample_b0(double *b) {
  double startlike;
  if ( verbose>1 ) {
    startlike = likelihood();
    yap_message("sample_b0 (pre): b0=%lf, lp=%lf\n",
                *b, startlike);
  }
  if ( ddP.PYbeta==H_HDP )
    sample_b0_DP(&ddP.b0);
  else 
    sample_b0_PDD(&ddP.b0);
#if 0
  if ( *b<PYP_CONC_MIN )
    *b = PYP_CONC_MIN;
  if ( *b>PYP_CONC_MAX )
    *b = PYP_CONC_MAX;
#endif
  cache_update("b0");
  if ( verbose>1 ) {
    double endlike = likelihood();
    yap_message("sample_b0 (post): b0=%lf, lp=%lf\n",
                *b, endlike);
    if ( 0 && endlike < startlike-50 )
      yap_quit("Sampler failed due to huge decrease: from %lf to %lf!\n",
	       startlike, endlike);
  }
}
