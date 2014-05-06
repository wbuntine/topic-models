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

//  #define A_DEBUG

#ifdef A_DEBUG
  static double last_val = 0;
  static double last_like = 0;
#endif
/*
 *  globals from hca.h and stats.h
 */
extern int verbose;
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


