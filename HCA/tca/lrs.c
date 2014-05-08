/*
 * Various likelihood calculations
 * Copyright (C) 2011-2014 Wray Buntine
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@monash.edu)
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "yap.h"
#include "util.h"
#include "stable.h"
#include "lgamma.h"
#include "tca.h"
#include "data.h"
#include "stats.h"


/*
 *   does two types of likelihood calcs:
 *
 *     GibbsNone:  adds own counts to the populations to
 *                 affect the estimated model ... but only
                   temporarily during sampling
 *     GibbsHold:  splits data into 2 ... first part is like
 *                 GibbsNone, second part uses the topic proportions
 *                 estimated from first and then does a standard
 *                 word-prob. estimate using the topic proportions
 *
 *   NB.  topic assignments z[] left as they are, so if
 *        called again have been warmed up
 */
double lp_test_ML(/*
		   *   fix==GibbsHold for hold-out testing
		   *   fix==GibbsNone for old max. likelihood testing
		   */
		  enum GibbsType fix)
{
  int i, r;
  float *fact = fvec(ddN.T);
  int StartTestDoc=ddN.D-ddN.TEST, EndTestDoc=ddN.D;
  double lik=0.0;
  int totw=0;
  D_MiSi_t dD;

  /*
   *  must account for other totals over docs
   *  which would be modified by adding the test doc
   *    m_evt[], Nt[] - we don't change these, usd in wordfact()
   *    C_eDt[t]  = \sum_d c_dt[d][t]
   *    C_e = \sum_t C_eDt[t]
   *
   *   assume topic assignments are set up in z[], 
   *   but stats not added elsewhere
   */
  // check_m_evt(0);
  if ( PCTL_BURSTY() ) misi_init(&ddM,&dD);

  /*
   *   now run sampler on all test docs
   */
  for(i=StartTestDoc; i<EndTestDoc; i++) {
    double hmean = -1e30;
#ifdef TRACE_WT
    int e = ddD.e[i];
#endif
    int  thisw =  add_doc(i, fix);
    if ( thisw<=1 || (fix==GibbsHold && thisw>=ddD.N_dT[i]-1) ) {
#ifdef TRACE_WT
      yap_message("remove_doc(d=%d,N=%d,T=%d) before continue\n",
		  i, (int)ddS.m_evt[e][TR_W][TR_T],(int)ddS.s_evt[e][TR_W][TR_T]);
#endif
      remove_doc(i, fix);
#ifdef TRACE_WT
      yap_message("after remove_doc(d=%d,N=%d,T=%d)\n",
		  i, (int)ddS.m_evt[e][TR_W][TR_T],(int)ddS.s_evt[e][TR_W][TR_T]);
#endif
      continue;
    }
    if ( PCTL_BURSTY() ) misi_build(&dD,i,0);
    for (r=0; r<ddP.mltburn; r++) 
      gibbs_lda(fix, i, ddD.N_dT[i], fact, &dD);
    /*
     *   record harmonic mean of last (samples-burnin)
     */
    for (; r<ddP.mltiter; r++) 
      hmean = logadd(hmean,-gibbs_lda(fix, i, ddD.N_dT[i], fact, &dD));
    lik += log(ddP.mltiter-ddP.mltburn) - hmean;
    totw += thisw;
    if ( PCTL_BURSTY() ) misi_unbuild(&dD,i,0);
#ifdef TRACE_WT
    yap_message("remove_doc(d=%d,N=%d,T=%d) end loop\n",
		i, (int)ddS.m_evt[e][TR_W][TR_T],(int)ddS.s_evt[e][TR_W][TR_T]);
#endif
    remove_doc(i, fix);
#ifdef TRACE_WT
    yap_message("after remove_doc(d=%d,N=%d,T=%d) end loop\n",
		i, (int)ddS.m_evt[e][TR_W][TR_T],(int)ddS.s_evt[e][TR_W][TR_T]);
#endif
    //  yap_message("%d:  %lf-%lf / %d\n", i, hmean, log(ddP.mltiter-ddP.mltburn) - hmean, thisw); 
    // check_m_evt(ddD.e[i]);
  }
  free(fact);
  if ( PCTL_BURSTY() ) misi_free(&dD);
  if ( totw==0 )
    return 0;
  return lik/totw;
}



