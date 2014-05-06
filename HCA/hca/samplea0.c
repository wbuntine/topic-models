/*
 * Sampling utility for a0
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
#include "stats.h"

//  #define A0_DEBUG

#ifdef A0_DEBUG
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
#ifdef A0_DEBUG
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
#ifdef A0_DEBUG
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

void sample_a0(double *mya0) {
#ifdef A0_DEBUG
  last_val = 0;
  last_like = 0;
#endif  
  if ( verbose>1 )
    yap_message("sample_a0 (pre):  a0=%lf, lp=%lf\n",
		*mya0, likelihood());
  myarms(PYP_DISC_MIN, PYP_DISC_MAX, &a0terms, NULL, mya0, "a0");
  cache_update("a0");
  if ( verbose>1 )
    yap_message("sample_a0 (post):  a0=%lf, lp=%lf\n",
		*mya0, likelihood());
}


