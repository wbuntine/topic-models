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

/*
 *     these are set in the sample_b() routine
 *     for use in bterm() 
 */
static long Tsum;
static uint16_t *localTd;

/*
 *
 */
static double bterms(double b, void *mydata) {
  int i;
  double val = pctl_gammaprior(b);
  double lgb = lgamma(b);
  double lgba;
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
 *    this is the sampler given in Lan Du's papers
 */
void sample_b(double *b) {
  double startlike;
  int i, t;
  if ( verbose>1 ) {
    startlike = likelihood();
    yap_message("sample_b (pre): b=%lf, lp=%lf\n",
		*b, startlike);
  }
  /*
   *   compute t totals for docs since not stored
   */
  localTd = malloc(sizeof(*localTd)*ddN.DT);
  for (i=0; i<ddN.DT; i++) {
    uint16_t Td_ = 0;
    for (t=0; t<ddN.T; t++)
      Td_ += ddS.Tdt[i][t];
    Tsum += Td_;
    localTd[i] = Td_;
  }

  myarmsMH(PYP_CONC_MIN, PYP_CONC_MAX, &bterms, NULL, b, "b", 1);
  cache_update("b");
  if ( verbose>1 ) {
    double endlike = likelihood();
    yap_message("sample_b (post): b=%lf, lp=%lf\n",
		*b, endlike);
    if ( endlike < startlike-50 ) {
      yap_quit("Sampler failed due to huge decrease!\n");
    }
  }
  free(localTd);
}

 
