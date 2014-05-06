/*
 * Sampling utility for bw
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
 *   NOTE:    the auxillary variable method introducing
 *            "q" with a beta distribution means that all b>0,
 *            although a PDP/PDD allows b>-a
 *            i.e.,   our method restricts "b" a little
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


// #define BW_DEBUG

// #define SBW_USECACHE
/*
 *
 */
static double bwterms(double bw, void *mydata) {
  int t;
  double val = pctl_gammaprior(bw);
#ifdef SBW_USECACHE
  struct gcache_s lgba_t;
  struct gcache_s lgb_t;
#else
  double lgb = lgamma(bw);
  double lgba;
#endif
#ifdef SBW_USECACHE
  if ( ddP.awpar>0 )
    gcache_init(&lgba_t, bw/ddP.awpar);
  gcache_init(&lgb_t, bw);
#else
  if ( ddP.awpar>0 )
    lgba = lgamma(bw/ddP.awpar);
#endif
  for (t=0; t<ddN.T; t++) {
#ifdef SBW_USECACHE
    val += gcache_value(&lgba_t, (int)ddS.TWt[t]) 
      - gcache_value(&lgb_t, (int)ddS.NWt[t]);
#else
    if ( ddP.awpar>0 )
      val += gammadiff(ddS.TWt[t], bw/ddP.awpar, lgba);
    else 
      val += ddS.TWt[t] * log(bw);
    val -= gammadiff((int)ddS.NWt[t], bw, lgb) ;
#endif
  }
  myarms_evals++;
#ifdef BW_DEBUG
  yap_message("Eval bwterms(%lf) = %lf", bw, val);
  ddP.bwpar = bw;
  cache_update("bw");
  yap_message(", lp=%lf\n", likelihood());
#endif
  return val;
}

/*
 *    this is the sampler given in Lan Du's papers
 */
void sample_bw(double *bw) {
  double startlike;
  if ( verbose>1 ) {
    startlike = likelihood();
    yap_message("sample_bw (pre): bw=%lf, lp=%lf\n",
		*bw, startlike);
  }
  myarmsMH(PYP_CONC_MIN, PYP_CONC_MAX, &bwterms, NULL, bw, "bw", 1);
  cache_update("bw");
  if ( verbose>1 ) {
    double endlike = likelihood();
    yap_message("sample_bw (post): bw=%lf, lp=%lf\n",
		*bw, endlike);
    if ( endlike < startlike-50 ) {
      int i;
      yap_message("K = %d;\n", ddN.T);
      yap_message("T[%d] = {", ddN.T);
      for (i=0; i<ddN.T; i++) {
	yap_message("%d", (int)ddS.TWt[i]);
	if ( i<ddN.T-1 )
	  yap_message(",");
      }
      yap_message("};\nN[%d] = {", ddN.T);
      for (i=0; i<ddN.T; i++) {
	yap_message("%d", (int)ddS.NWt[i]);
	if ( i<ddN.T-1 )
	  yap_message(",");
      }
      yap_message("};\naw=%lf;\n", ddP.awpar);
      yap_message("bw=%lf;\n", ddP.bwpar);
      yap_quit("Sampler failed due to huge decrease!\n");
    }
  }
}

