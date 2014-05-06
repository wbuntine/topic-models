/*
 * Sampling utility for alpha
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

#include "util.h"
#include "yap.h"
#include "lgamma.h"
#include "myarms.h"
#include "sample.h"
#include "hca.h"
#include "stats.h"

#define CONJPRIOR
//  #define ALPHA_DEBUG

#ifdef ALPHA_DEBUG
  static double last_val = 0;
  static double last_like = 0;
#endif
/*
 */
static double alphaterms(double alpha, void *mydata) {
  int t,s;
  double val = 0;
  double tot;
  double lga = lgamma(alpha);
  double lgat = lgamma(ddN.T*alpha);
#ifdef ALPHA_DEBUG
  double like;
#endif
#ifdef CONJPRIOR
  val += ddN.T*(lgamma(alpha+1.0/ddN.T) - lga);
  val -= lgamma(ddN.T*alpha+1) - lgamma(ddN.T*alpha);
#endif
  for (s=0; s<ddN.DT; s++) {
    tot = 0;
    for (t=0; t<ddN.T; t++) {
      tot += alpha+ddS.Ndt[s][t];
      val += gammadiff(ddS.Ndt[s][t],alpha,lga);
    }
    val -= lgamma(tot) - lgat;
  }
  myarms_evals++;
  myarms_last = alpha;
#ifdef ALPHA_DEBUG
  yap_message("Eval alphaterms(%lf) = %lf\n", alpha, val);
  ddP.alpha = alpha;
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
 *  assumes uniform prior Dirichlet
 */
void sample_alpha(double *alpha) {
#ifdef ALPHA_DEBUG
  last_val = 0;
  last_like = 0;
#endif  
  if ( verbose>1 ) 
    yap_message("sample_alpha (pre):  alpha=%lf, lp=%lf\n",
		*alpha, likelihood());
  if ( myarmsMH(DIR_MIN, DIR_MAX, 
		&alphaterms, NULL, alpha, "alpha",1) ) {
    yap_message("sample_alpha: error in result\n");
  } 
  cache_update("alpha");
  if ( verbose>1 ) 
    yap_message("sample_alpha (post):  alpha=%lf, lp=%lf\n",
		*alpha, likelihood());
}

