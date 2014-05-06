/*
 * Lower level routines for changing stats
 * Copyright (C) 2010-2013 Wray Buntine 
 *           (C) 2014 Swapnil Mishra and Wray Buntine
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@monash.edu)
 *         Swapnil Mishra (swapnil.mishra@anu.edu.au)
 *
 *     Various data structure read/write/report routines.
 *     Defined in "hca.h"
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "yap.h"
#include "util.h"
#include "hca.h"
#include "data.h"
#include "stats.h"
#include "pctl.h"
#include "atomic.h"

/*
 *   basically, we abandon all stats for this doc
 *   better have not had other totals in ddS.TDt, ddS.TWTmz, etc.
 */
void zero_doc(int d) {
  int t;
  for (t=0; t<ddN.T; t++) {
    ddS.Ndt[d][t] = 0;
  }
  ddS.NdT[d] = 0;
  if ( ddP.PYalpha ) {
    for (t=0; t<ddN.T; t++) {
      ddS.Tdt[d][t] = 0;
    }
  }
}

void unfix_tableidtopic(int d, int t) {
  int val;
  ddS.Tdt[d][t]--;
  assert(ddS.Tdt[d][t]>=0);
  assert(ddS.Tdt[d][t]>0 || ddS.Ndt[d][t]==0);
  val = atomic_decr(ddS.TDt[t]);
  atomic_decr(ddS.TDT);
  if ( val==0 ) {
    atomic_decr(ddS.TDTnz);
    ddS.Tlife[t] = 0;
  }
#ifdef CACHE_ABTP
  alphabasetopicprob(-(ddN.T+1));
#endif
}

void fix_tableidtopic(int d, int t) {
  int val;
  assert(ddS.Tdt);
  ddS.Tdt[d][t]++;
  val = atomic_incr(ddS.TDt[t]);
  if ( val==1 )
    atomic_incr(ddS.TDTnz);
  atomic_incr(ddS.TDT);
#ifdef CACHE_ABTP
  alphabasetopicprob(-(ddN.T+1));
#endif
}

void unfix_tableidword(int w, int t) {
  int val;
  assert(ddS.Twt[w][t]);
  atomic_decr(ddS.Twt[w][t]);
  val = atomic_decr(ddS.TwT[w]);
  atomic_decr(ddS.TWt[t]);
  atomic_decr(ddS.TWT);
  if ( val==0 ) {
    atomic_decr(ddS.TWTnz);
  }
}

void fix_tableidword(int w, int t) {
  int val;
  atomic_incr(ddS.Twt[w][t]);
  atomic_incr(ddS.TWt[t]);
  val = atomic_incr(ddS.TwT[w]);
  if ( val==1 ) {
    atomic_incr( ddS.TWTnz);
  }
  atomic_incr(ddS.TWT);
}

/*
 *    remove affects of document from stats
 */
int remove_doc(int d, enum GibbsType fix) {
  int i, t;
  for (t=0; t<ddN.T; t++) 
    ddS.Ndt[d][t] = 0;
  ddS.NdT[d] = 0;
  if ( ddP.PYalpha ) {
#ifdef CACHE_ABTP
    alphabasetopicprob(-(ddN.T+1));
#endif
    for (t=0; t<ddN.T; t++) 
      if ( ddS.Tdt[d][t]>0 ) {
        int val;
	val = atomic_sub(ddS.TDt[t],ddS.Tdt[d][t]);
	atomic_sub(ddS.TDT,ddS.Tdt[d][t]);
	if ( val==0 ) {
          atomic_decr(ddS.TDTnz);
          ddS.Tlife[t] = 0;
        }
	ddS.Tdt[d][t] = 0;
      }
  }
  for (i=ddD.NdTcum[d]; i<ddD.NdTcum[d+1]; i++) {
    if ( fix!=GibbsHold || !pctl_hold(i) ) {
      /*
       *   these words are for training
       */
      if ( (ddP.bdk==NULL) || Z_issetr(ddS.z[i]) ) {
	t = Z_t(ddS.z[i]);
	if ( ddP.phi==NULL ) {
          int val;
	  int w = ddD.w[i];
	  atomic_decr(ddS.NWt[t]);
	  // assert(ddS.Nwt[w][t]>0);
	  val = atomic_decr(ddS.Nwt[w][t]);
	  if (  ddP.PYbeta ) 
	    if ( val==0 || ddS.Twt[w][t]>ddS.Nwt[w][t] ) 
	      unfix_tableidword(w,t);  // ???? WRONG ??
	}
      }
    }
  }
  return 0;
}

/*
 *    add affects of document to stats
 *
 *    when using multis, reset ddM.Mi[] to be totals
 *    for training words only, ignore test words in hold
 */
int add_doc(int d, enum GibbsType fix) {
  int i, t, w, nd=0;
  int mi = 0;
  if ( ddP.bdk!=NULL )
    mi = ddM.MI[d];

  for (i=ddD.NdTcum[d]; i<ddD.NdTcum[d+1]; i++) {
    if ( fix!=GibbsHold || pctl_hold(i) ) {
       /*
       *   these words are for perp. calcs
       */     
      nd++;
    }
    if ( fix!=GibbsHold || !pctl_hold(i) ) {
      t = Z_t(ddS.z[i]);
      /*
       *   these words are for training
       */
      ddS.Ndt[d][t]++;
      ddS.NdT[d]++;	
      if ( (ddP.bdk==NULL) || Z_issetr(ddS.z[i]) ) {
	if ( ddP.phi==NULL ) {
          int val;
	  w = ddD.w[i];
	  atomic_incr(ddS.NWt[t]);
	  val = atomic_incr(ddS.Nwt[w][t]);
	  if (  ddP.PYbeta && val==1 ) 
	    fix_tableidword(w,t);
	}
      }
    }
    if ( (ddP.bdk!=NULL) && M_multi(i) ) mi++;
  }
  if ( ddP.PYalpha) {
    /*  initialise ddS.Tdt[d][*]  */ 
#ifdef CACHE_ABTP
    alphabasetopicprob(-(ddN.T+1));
#endif
    /*
     *   adjust table count stats based on Ndt[d]
     */
    for (t=0; t<ddN.T; t++) {
      if ( ddS.Ndt[d][t]>0 ) {
        int val;
	ddS.Tdt[d][t] = 1;
        val = atomic_incr(ddS.TDt[t]);
        if ( val==1 )
          atomic_incr(ddS.TDTnz);
        atomic_incr(ddS.TDT);
      } else
	ddS.Tdt[d][t] = 0;
    }
  }
  // yap_message("add_doc(%d):  %d\n", d, nd);
  return nd;
}

/*
 *  coount non-zero entries in Nwt[][t]
 */
int nonzero_Nwt(int t) {
  int w;
  int nz=0;
  for (w=0; w<ddN.W; w++ )
    if ( ddS.Nwt[w][t]>0 )
      nz++;
  return nz;
}

/*
 *  count non-zero entries in Ndt[][t]
 */
int nonzero_Ndt(int t) {
  int d;
  int nz=0;
  for (d=0; d<ddN.DT; d++ )
    if ( ddS.Ndt[d][t]>0 )
      nz++;
  return nz;
}

/*
 *   find total t's in document
 */
uint16_t comp_Td(int did) {
  uint16_t Td_ = 0;
  int t;
  assert(ddS.Tdt);
  for (t=0; t<ddN.T; t++) {
    Td_ += ddS.Tdt[did][t];
  }
  return Td_;
}
