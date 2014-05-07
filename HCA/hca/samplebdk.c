/*
 * Sampling utility for bdk
 * Copyright (C) 2013 Wray Buntine 
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

/*
 *     define to make sample independently for different k
 */
#define SAMPLEBDK_MULTI

/*
 *    docstats[d] points to store for doc d
 *    
 *    stores blocks per nonzero topic
 *       [0] = topic
 *       [1] = total count
 *       [2] = total tables
 *
 *    ending block
 *       [0] = ddN.T+1
 */
static uint16_t **docstats;

/*
 *  globals from hca.h and stats.h
 */
extern double likelihood_bdk();

// #define B_DEBUG

/*
 *    just call the likelihood function
 */
#ifdef SAMPLEBDK_MULTI
static double bdkterms(double b, void *mydata) {
  int t = *((int *)mydata);
  double val;
  ddP.bdk[t] = b;
  val = dmi_likelihood_bterms(&ddM, t, docstats, pctl_gammaprior,
			      ddP.ad, ddP.bdk);
  myarms_evals++;
#ifdef B_DEBUG
  yap_message("Eval (from likelihood) bdk[%d]terms(%lf) = %lf", t,b,val);
#endif
  return val;
}
#else
static double bdkterms(double b, void *mydata) {
  int t;
  double val;
  for (t=0; t<ddN.T; t++) 
    ddP.bdk[t] = b;
  val = dmi_likelihood_bterms(&ddM, -1, docstats, pctl_gammaprior,
			      ddP.ad, ddP.bdk);
  myarms_evals++;
#ifdef B_DEBUG
  yap_message("Eval (from likelihood) bdkterms(%lf) = %lf\n", b, val);
#endif
  return val;
}
#endif

/*
 *    this uses ARMS
 */
void sample_bdk(double *b) {
#ifdef SAMPLEBDK_MULTI
  static int batchstart = 0;
  int tt;
#endif
  double startlike;
  int t;
  assert(b);
  docstats = dmi_bstore(&ddM);
#ifdef SAMPLEBDK_MULTI
  assert(ddP.bdkbatch>0);
  if ( verbose>1 ) {
      t = batchstart%ddN.T;
      startlike = likelihood();
      yap_message("sample_bdk[%d:%d] (pre): b=%lf, lp=%lf\n",
		  t, (ddP.bdkbatch+batchstart-1)%ddN.T,b[t], startlike);
  }
  for (tt=0; tt<ddP.bdkbatch; tt++) {
    t = (batchstart+tt)%ddN.T;
    
    myarmsMH(PYP_CONC_MIN, PYP_CONC_MAX, &bdkterms, &t, &b[t], "bdk", 1);
  }
  if ( verbose>1 ) {
      double endlike = likelihood();
      t = (batchstart)%ddN.T;
      yap_message("sample_bdk[%d:%d] (post): bdk=%lf, lp=%lf\n",
		  t, (ddP.bdkbatch+batchstart-1)%ddN.T, b[t], endlike);
      if ( endlike < startlike-50 ) {
	yap_quit("Sampler failed due to huge decrease!\n");
      }
  }
  batchstart = (batchstart+tt)%ddN.T;
#else    
  if ( verbose>1 ) {
    startlike = likelihood();
    yap_message("sample_bdk (pre): bdk[*]=%lf, lp=%lf\n",
		b[0], startlike);
  }
    
  myarmsMH(PYP_CONC_MIN, PYP_CONC_MAX, &bdkterms, NULL, &b[0], "bdk", 1);
  for (t=1; t<ddN.T; t++)
    b[t] = b[0];
  if ( verbose>1 ) {
    double endlike = likelihood();
    yap_message("sample_bdk (post): bdk[*]=%lf, lp=%lf\n",
		b[0], endlike);
    if ( endlike < startlike-50 ) {
      yap_quit("Sampler failed due to huge decrease!\n");
    }
  }
#endif
  dmi_freebstore(&ddM,docstats);
}

 
