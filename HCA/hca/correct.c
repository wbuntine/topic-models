/*
 * Various data structures for statistics
 * Copyright (C) 2010-2014 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@monash.edu)
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

/*
 *  ensure Twt satisfies constraints
 */
void correct_twt()  {
  int w, t;
  
  if ( ddP.PYbeta==0 )
    return;

  ddS.TWT = ddS.TWTnz = 0;
  memset((void*)ddS.TwT, 0, sizeof(ddS.TwT[0])*ddN.W);
  memset((void*)ddS.TWt, 0, sizeof(ddS.TWt[0])*ddN.T);
  for(t = 0; t < ddN.T; t++) {
    for(w = 0; w < ddN.W; w++) {
      if ( ddS.Twt[w][t]>ddS.Nwt[w][t] )
	ddS.Twt[w][t] = ddS.Nwt[w][t];
      if ( ddS.Twt[w][t]==0 && ddS.Nwt[w][t]>0 )
	ddS.Twt[w][t] = 1;
      ddS.TWt[t] += ddS.Twt[w][t];
      ddS.TwT[w] += ddS.Twt[w][t];
    }
    ddS.TWT += ddS.TWt[t];
  }
  for(w = 0; w < ddN.W; w++) {
    if ( ddS.TwT[w]>0 )
      ddS.TWTnz ++;
  }
}

void correct_tdt(int reset)  {
  int d, t;
  
  if ( ddP.PYalpha==0 )
    return;

#ifdef DEBUGTEST
  if ( reset==0 ) {
    uint16_t **Tdt = u16mat(ddN.D,ddN.T);
    int i;
    for (d = 0; d < ddN.DT; d++) {
      for (i=ddD.NdTcum[d]; i<ddD.NdTcum[d+1]; i++) 
	Tdt[d][Z_t(ddS.z[i])]++;
      for (t = 0; t < ddN.T; t++) {
	if ( ddS.Tdt[d][t]!=Tdt[d][t] )
	  yap_message("Unequal Tdt totals for doc %d\n", d);
      }
    }
    free(Tdt[0]); free(Tdt);
  }
#endif

  /*
   *  reset Tdt
   */
  if ( reset ) {
    for (d=0; d<ddN.D; d++) {
      int i;
      /*   ddS.Tdt not allocated monolithically  */
      memset((void*)ddS.Tdt[d], 0, sizeof(ddS.Tdt[0][0])*ddN.T);
      for (i=ddD.NdTcum[d]; i<ddD.NdTcum[d+1]; i++) 
	ddS.Tdt[d][Z_t(ddS.z[i])]++;
    }
  }
  /*
   *  correct derived stats
   */
  ddS.TDT = 0;
  ddS.TDTnz = 0;
  memset((void*)ddS.TDt, 0, sizeof(ddS.TDt[0])*ddN.T);
  for (t = 0; t < ddN.T; t++) {
    for (d = 0; d < ddN.D; d++) {
      if ( reset==0 ) {
	if ( ddS.Tdt[d][t]>ddS.Ndt[d][t] )
	  ddS.Tdt[d][t] = ddS.Ndt[d][t];
	if ( ddS.Tdt[d][t]==0 && ddS.Ndt[d][t]>0 )
	  ddS.Tdt[d][t] = 1;
      }
      if ( d<ddN.DT ) 
	ddS.TDt[t] += ddS.Tdt[d][t];
    }
    ddS.TDT += ddS.TDt[t];
  }
  for(t = 0; t < ddN.T; t++) {
    if ( ddS.TDt[t]>0 )
      ddS.TDTnz ++;
  }
}

#ifdef NG_SCALESTATS
/*
 *   occasionally recompute them
 */
void correct_NGscalestats(int redo) {
  if ( redo || ++ddS.NGscalestats_recomp>1000 ) {
    int i, t;
    for (t=0; t<ddN.T; t++) 
      ddS.NGscalestats[t] = 0;
    for (i=0; i<ddN.DT; i++) {
      if ( ddS.NdT[i]==0 || ddS.UN[i]==0 ) continue;
      for (t=0; t<ddN.T; t++) {
	if ( !ddS.sparse || M_docsparse(i,t) )
	  ddS.NGscalestats[t] += log(1.0+ddS.UN[i]/ddP.NGbeta[t]);
      }
    }
    ddS.NGscalestats_recomp = 0;
  }
}
#endif

void correct_docsp() {
  int i, t;
  if ( ddS.sparse==NULL )
    return;
  /*
   *  ensure non-zero topics are not sparse
   */
  for (i=0; i<ddN.D; i++) {
    if ( ddS.NdT[i]==0 || ddS.UN[i]==0 )
      continue;
    for (t=0; t<ddN.T; t++) 
      if ( ddS.Ndt[i][t]>0 )
	M_docsp_set(i,t);
  }
  /*
   *  fix up total vectors
   */
  for (t=0; t<ddN.T; t++) {
    ddS.sparseD[t] = 0;
  }
  for (i=0; i<ddN.DT; i++) {
    if ( ddS.NdT[i]==0 || ddS.UN[i]==0 )
      continue;
    for (t=0; t<ddN.T; t++) {
      if ( M_docsparse(i,t) )
	ddS.sparseD[t]++;
    }
  }
#ifdef NG_SCALESTATS
  correct_NGscalestats(1);
#endif
}
