/*
 * Sampling utility for beta
 * Copyright (C) 2010 Wray Buntine 
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
#include "util.h"
#include "lgamma.h"
#include "myarms.h"
#include "sample.h"
#include "stats.h"
#include "pctl.h"

//  #define BETA_DEBUG

#ifdef BETA_DEBUG
  static double last_val = 0;
  static double last_like = 0;
#endif

/*
 *  globals from lrs.
 */
extern int verbose;

/*
 */
static double old_beta;

/*
 */
static double betaterms(double mytbeta, void *mydata) {
  int j,t;
  double val = 0;    
#ifdef BETA_DEBUG
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
#ifdef BETA_DEBUG
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

/*
 *  this allows Dirichlet prior to be non-uniform,
 *  so optimisation done on total weight;
 */
void sample_beta(double *mytbeta) {
  double bmax =  DIR_MAX*ddN.W;
#ifdef BETA_DEBUG
  last_val = 0;
  last_like = 0;
#endif  
  // if ( bmax > DIR_TOTAL_MAX )  bmax = DIR_TOTAL_MAX;
  /*
   *  beta is set using ddP.beta, so need to factor it out
   */
  old_beta = ddP.beta;
  if ( verbose>1 )
    yap_message("sample_beta (pre):  beta=%lf, lp=%lf\n",
		*mytbeta, likelihood());
  myarmsMH(DIR_MIN*ddN.W, bmax, &betaterms, NULL, mytbeta, "beta",1);
  cache_update("beta");
  if ( verbose>1 ) {
    yap_message("sample_beta (post): beta=%lf, lp=%lf\n",
		*mytbeta, likelihood());
  }
}
