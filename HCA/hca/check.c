/*
 * Checking routines - debugging only
 * Copyright (C) 2010-2013 Wray Buntine 
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

void check_Tw() {
#ifndef NDEBUG
  int w;
  int ttot=0, nztot=0, ztot=0;
  for (w=0; w<ddN.W; w++) {
    if ( ddS.TwT[w]>0 )
      nztot++;
    else
      ztot++;
    ttot += ddS.TwT[w];
  }
  if ( ttot != ddS.TWT ) 
    yap_message(" \\sum_w TwT[w] != TWT\n"); 
  if ( nztot != ddS.TWTnz ) {
    yap_message(" \\sum_w 1(TwT[w]) != TWTnz\n"); 
    yap_message("check_Tw() z=%d, nz=%d, W=%d\n",
		ztot, nztot, ddN.W);  
  }
#endif
}

void check_Nwt(int w, int val) {
#ifndef NDEBUG
  int t;
  int ntot=0, ttot=0;
  if ( w>=ddN.W )
    return;
  for (t=0; t<ddN.T; t++) {
    ntot += ddS.Nwt[w][t];
    if ( ddS.Twt )
      ttot += ddS.Twt[w][t];
  }
  if ( val>0 && ntot!=val ) {
    yap_message(" (NwT[%d]=%d,NwT[%d]=%d)", w, ntot, w, ttot);
  }
#endif
}

void check_Ndt(int d) {
#ifndef NDEBUG
   int t;
  for (t=0; t<ddN.T; t++) {
    if ( ddS.Ndt[d][t]==0 && ddS.Tdt[d][t]>0 )
      assert(ddS.Ndt[d][t]>0 || ddS.Tdt[d][t]==0);
    if ( ddS.Ndt[d][t]>0 && ddS.Tdt[d][t]==0 )
      assert(ddS.Ndt[d][t]==0 || ddS.Tdt[d][t]>0);
  }
#endif
}

#ifdef NG_SPARSE
void check_sparse() {
#ifndef NDEBUG
  int i,k;
  int cnt;
  for (k=0; k<=ddN.T; k++) {
    cnt = 0;
    for (i=0; i<ddN.DT; i++) {
      if ( M_docsparse(i,k) ) {
	cnt++;
	if ( ddS.Ndt[i][k]==0 )
	    yap_quit("ddS.sparse[%d][%d]==1 but actual count 0 in check_sparse()\n", i, k);
      }
    }
    if ( cnt!=ddS.sparseD[k] )
      yap_quit("ddS.sparseD[%d]!=actual in check_sparse()\n", k);
  }
#endif
}
#endif

void check_TWT() {
#ifndef NDEBUG
  int w;
  int tot = 0;
  int totnz = 0;
  for (w=0; w<ddN.W; w++) {
    tot += ddS.TwT[w];
    if ( ddS.TwT[w]>0 ) totnz++;
  }
  if ( tot != ddS.TWT )
    yap_message("ddS.TWT = %d, computed %d\n", ddS.TWT, tot);
  if ( totnz != ddS.TWTnz )
    yap_message("ddS.TWTnz = %d, computed %d\n", ddS.TWTnz, totnz);
#endif
}
