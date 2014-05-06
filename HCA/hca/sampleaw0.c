/*
 * Sampling utility for aw0
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

//  #define AW0_DEBUG

#ifdef AW0_DEBUG
  static double last_val = 0;
  static double last_like = 0;
#endif
/*
 *  globals from hca.h and stats.h
 */
extern int verbose;

/*
 */
static double aw0terms(double myaw0, void *mydata) {
  int i;
  double l1aw0 = log(1-myaw0);
  double l2aw0 = log((1-myaw0)*(2-myaw0));
  double lgaw0 = lgamma(1-myaw0);
  double val = 0;
#ifdef AW0_DEBUG
  double like;
#endif
  val += ddS.TWTnz*log(myaw0) + lgamma(ddP.bw0/myaw0+ddS.TWTnz) 
    - lgamma(ddP.bw0/myaw0);
  for (i=0; i<ddN.W; i++) 
    /*  note the root node is a PDD so all t's = 1 */
    if ( ddS.TwT[i]>1 ) {
      if ( ddS.TwT[i]==2 )
	val += l1aw0;
      else if ( ddS.TwT[i]==3 )
	val += l2aw0;
      else
	val += lgamma(ddS.TwT[i]-myaw0) - lgaw0;
    }
#ifdef AW0_DEBUG
  yap_message("Eval aw0terms(%lf) = %lf", myaw0, val);
  ddP.aw0 = myaw0;
  cache_update("aw0");
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

void sample_aw0(double *myaw0) {
#ifdef AW0_DEBUG
  last_val = 0;
  last_like = 0;
#endif  
  /*
   *   compute it in first pass,
   *   then use it inside aw0terms() and aw0terms_da()
   */
  if ( verbose>1 )
    yap_message("sample_aw0 (pre):  aw0=%lf, lp=%lf\n",
		*myaw0, likelihood());
  myarms(PYP_DISC_MIN, PYP_DISC_MAX, &aw0terms, NULL, myaw0, "aw0");
  cache_update("aw0");
  if ( verbose>1 )
    yap_message("sample_aw0 (post):  aw0=%lf, lp=%lf\n",
		*myaw0, likelihood());
}


