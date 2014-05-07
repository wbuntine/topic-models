/*
 * Sampling utility for adk
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

//  #define A_DEBUG

#ifdef A_DEBUG
  static double last_val = 0;
  static double last_like = 0;
#endif
/*
 *  globals from hca.h and stats.h
 */
extern int verbose;
double likelihood_bdk();
  
/*
 *    docstats[d] points to store for doc d
 *    
 *    stores blocks per nonzero topic
 *       [0] = topic
 *       [1] = total count
 *       [2] = total tables
 *       [3] = no. of Mi[] Si[] pairs following with Mi[]>=2
 *       [4+2w]+[5+2w] = Mi[]+Si[] for w-th word
 *
 *    ending block
 *       [0] = ddN.T+1
 */
static uint16_t **docstats;

/*
 */
static double adkterms(double mya, void *mydata) {
  double val = 0;
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

void sample_adk(double *mya) {
#ifdef A_DEBUG
  last_val = 0;
  last_like = 0;
#endif
  docstats = dmi_astore(&ddM);
  if ( verbose>1 )
    yap_message("sample_adk (pre):  ad=%lf, lp=%lf\n",
		*mya, likelihood());
  myarms(PYP_DISC_MIN, PYP_DISC_MAX, &adkterms, NULL, mya, "adk");
  cache_update("ad");
  dmi_freeastore(&ddM, docstats);
  docstats = NULL;
  if ( verbose>1 )
    yap_message("sample_adk (post):  ad=%lf, lp=%lf\n",
		*mya, likelihood());
 }

