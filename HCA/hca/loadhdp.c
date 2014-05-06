/*
 * Load the topic assignments file from Chang's HDP
 * Copyright (C) 2013 Wray Buntine 
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
#include <stdint.h>
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
 *    zero everything and rebuild entirely from Chang's HDP
 *            mode-word-assignments.dat
 *    but only for training docs
 *
 *    REALLY slow implementation
 */
#define HDP_FILE "mode-word-assignments.dat"
void hca_load_hdp(char *resstem) {  
  unsigned d, w, t, m;
  int last_d, last_m, last_t;
  int l;
  
  FILE *fp;
  
  /*
   *  initialisation *not* done for test docs
   */
  assert(ddP.PYbeta==H_None);
  assert(ddP.phi==NULL);

  for (l=0; l<ddN.NT; l++)
    ddS.z[l] = ddN.T+1;
  
  fp = popen("grep -v '^d w ' "  HDP_FILE " | sort -k1 -n -k3 -n -k4 -n -k3 -n ","r");
  if ( !fp ) 
    yap_quit("Cannot open '%s'\n", HDP_FILE);
  last_m = last_t = last_d = -1;
  while ( fscanf(fp, " %u %u %u %u", &d, &w, &t, &m)==4 ) {
    w--;
    if ( t>=ddN.T )
      yap_quit("Topic %d greater than set\n", t);
    if ( d>= ddN.DT ) 
      yap_quit("Word %d greater than set\n", w);
    if ( d != last_d && last_d>=0 ) {
      int tt;
      /*   
       *   new doc  
       */
      /*  check last one was filled properly */
      for (l=ddD.NdTcum[last_d]; l<ddD.NdTcum[last_d+1]; l++)
	if ( ddS.z[l]>=ddN.NT )
	  yap_quit("Document %d topic unset at %d-th position\n",
		   last_d, l);
      /*  check no docs skipped */
      if ( ddD.NdTcum[last_d+1] != ddD.NdTcum[d] )
	yap_quit("Document missed between %d and $d\n", last_d, d);
      /*  check Tdt[][] */
      for (tt=0; tt<ddN.T; tt++) {
	assert(ddS.Tdt[last_d][tt]<=ddS.Ndt[last_d][tt]);
      }
      last_m = last_t = -1;
    }
    for (l=ddD.NdTcum[d]; l<ddD.NdTcum[d+1]; l++) {
      if ( ddS.z[l]>=ddN.T && ddD.w[l]==w ) {
	/*
	 *  here is a matching word
	 */
	ddS.z[l] = t;
	ddS.NWt[t]++;
	ddS.Nwt[w][t]++;	
	ddS.Ndt[d][t]++;
	ddS.NdT[d]++;
	if ( last_d!=d ||
	     (last_d==d && last_t!=t) ||
	     (last_d==d && last_t==t && last_m!=m ) )
	  /*  new table for topic */
	  fix_tableidtopic(d,t);
	break;
      }
    }
    if ( l>=ddD.NdTcum[d+1] )
      yap_quit("No matching word %d found in doc %d \n", w+1, d);
    last_d = d;
    last_m = m;
    last_t = t;
  }
  for (t=0; t<ddN.T; t++) {
    // assert(ddS.Tdt[last_d][t]<=ddS.Ndt[last_d][t]);
  }
  for (l=0; l<ddN.NT; l++)
    if ( ddS.z[l]>= ddN.T ) 
      yap_quit("Document %d topic unset at %d-th position after finish\n",
	       ddD.d[l], l);
  pclose(fp);
}
