/*
 * Sampling utility for aw
 * Copyright (C) 2011-2103 Wray Buntine 
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

//  #define AW_DEBUG

#ifdef AW_DEBUG
  static double last_val = 0;
  static double last_like = 0;
#endif
/*
 *  globals from hca.h and stats.h
 */
extern int verbose;
/*
 */
static double awterms(double myaw, void *mydata) {
  int i, t;
  double val = 0;
  double law = log(myaw);
#ifdef AW_DEBUG
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
    val += Tw_*law + lgamma(ddP.bwpar/myaw+Tw_) - lgamma(ddP.bwpar/myaw);
  }  
  myarms_evals++;
#ifdef AW_DEBUG
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
#ifdef AW_DEBUG
  last_val = 0;
  last_like = 0;
#endif    /*
   *   compute it in first pass,
   *   then use it inside awterms() and awterms_da()
   */
  if ( verbose>1 )
    yap_message("sample_aw (pre):  aw=%lf, lp=%lf\n",
		*myaw, likelihood());
  myarms(PYP_DISC_MIN, PYP_DISC_MAX, &awterms, NULL, myaw, "aw");
  cache_update("aw");
  if ( verbose>1 )
    yap_message("sample_aw (post):  aw=%lf, lp=%lf\n",
		*myaw, likelihood());
 }


