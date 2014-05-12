/*
 * Deal with the caches and tables
 * Copyright (C) 2012-2014 Wray Buntine 
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
#include <string.h>
#include <math.h>
#include <assert.h>

#include "yap.h"
#include "yaps.h"
#include "lgamma.h"
#include "tca.h"
#include "data.h"
#include "pctl.h"

#define mymax(A,B) ((A>B)?(A):(B))

void cache_init() {
  /*
   *   now the libstb library has its own separate error handling,
   *   so make sure it goes through our yapper
   */
  yaps_yapper(yap_va);
 
  /*
   *   we'll share tables for all those with a==0,
   *   since it never changes
   */
  {
    char tag[200];
    int mN=0, mM=0, sN=0, sM=0;
    ddC.a_zero = NULL;
    tag[0] = 0;
    if ( ddP.a_mu==0 ) {
      mM = mymax(mM,ddP.maxM);
      mN = mymax(mN,ddP.maxN+1);
      sM = mymax(sM,100);
      sN = mymax(sN,ddP.maxN/5);
      strcat(tag, "a_mu, ");
    }
    if ( ddP.a_theta==0 ) {
      mM = mymax(mM,ddP.maxM);
      mN = mymax(mN,ddD.N_dTmax+1);
      sM = mymax(sM,100);
      sN = mymax(sN,ddD.N_dTmax/2);
      strcat(tag, "a_theta, ");
    }
    if ( ddP.a_phi0==0 ) {
      mM = mymax(mM,ddP.maxM);
      mN = mymax(mN,ddP.maxN+1);
      sM = mymax(sM,ddP.maxM/2);
      sN = mymax(sN,ddP.maxN/10);
      strcat(tag, "a_phi0, ");
    }
    if ( ddP.a_phi1==0 ) {
      mM = mymax(mM,ddP.maxM);
      mN = mymax(mN,ddP.maxN+1);
      sM = mymax(sM,ddP.maxM/2);
      sN = mymax(sN,ddP.maxN/10);
      strcat(tag, "a_phi1, ");
    }
    if ( PCTL_BURSTY() && ddP.a_burst==0 ) {
      mM = mymax(mM,ddP.maxM);
      mN = mymax(mN,ddD.N_dTmax+1);
      sM = mymax(sM,100);
      sN = mymax(sN,ddD.N_dTmax/2);
      strcat(tag, "a_burst, ");
    }
    if ( mM>0 ) {
      ddC.a_zero = S_make(sN, sM, mN, mM, 0,
#ifdef H_THREADS
                          S_THREADS|
#endif
                          S_STABLE|S_UVTABLE|S_FLOAT| S_QUITONBOUND);
      if ( ddP.a_mu==0 ) 
        ddC.a_mu = ddC.a_zero;
      if ( ddP.a_theta==0 )
        ddC.a_theta = ddC.a_zero;
      if ( ddP.a_phi0==0 ) 
        ddC.a_phi0 = ddC.a_zero;
      if ( ddP.a_phi1==0 ) 
        ddC.a_phi1 = ddC.a_zero;
      if ( PCTL_BURSTY() && ddP.a_burst==0 ) 
        ddC.a_burst = ddC.a_zero;
      strcat(tag, "all zero PYP");
      S_tag(ddC.a_zero, tag);
      if ( verbose ) 
        S_report(ddC.a_zero,NULL);
    }
  }
  if ( ddP.a_mu!=0 ) {
    ddC.a_mu = S_make(ddP.maxN/5, 100, ddP.maxN+1, ddP.maxM, ddP.a_mu, 
#ifdef H_THREADS
                      S_THREADS|
#endif
                      S_STABLE|S_UVTABLE|S_FLOAT| S_QUITONBOUND);
    S_tag(ddC.a_mu,"a_mu, topic time PYP");
    if ( verbose ) S_report(ddC.a_mu,NULL);
  }
  if ( ddP.a_theta!=0 ) {
    ddC.a_theta = S_make(ddD.N_dTmax/2, 100, ddD.N_dTmax+1, ddP.maxM, ddP.a_theta,
#ifdef H_THREADS
                         S_THREADS|
#endif
                         S_STABLE|S_UVTABLE|S_FLOAT| S_QUITONBOUND);
    S_tag(ddC.a_theta,"a_theta, topic->doc PYP");
    if ( verbose ) S_report(ddC.a_theta,NULL);
  }
  if ( ddP.a_phi0!=0 ) {
    ddC.a_phi0 = S_make(ddP.maxN/10, ddP.maxM/2, ddP.maxN+1, ddP.maxM, ddP.a_phi0, 
#ifdef H_THREADS
                        S_THREADS|
#endif
                        S_STABLE|S_UVTABLE|S_FLOAT| S_QUITONBOUND);
    S_tag(ddC.a_phi0,"a_phi0, word time PYP");
    if ( verbose ) S_report(ddC.a_phi0,NULL);
  }
  if ( ddP.a_phi1!=0 ) {
    ddC.a_phi1 = S_make(ddP.maxN/10, ddP.maxM/2, ddP.maxN+1, ddP.maxM, ddP.a_phi1, 
#ifdef H_THREADS
                        S_THREADS|
#endif
                        S_STABLE|S_UVTABLE|S_FLOAT| S_QUITONBOUND);
    S_tag(ddC.a_phi1,"a_phi1, word time PYP");
    if ( verbose ) S_report(ddC.a_phi1,NULL);
  }
  if ( !PCTL_BURSTY() ) {
    ddC.a_burst = NULL;
  } else if ( ddP.a_burst!=0 ) {
    ddC.a_burst = S_make(ddD.N_dTmax/2, 100, ddD.N_dTmax+1, ddP.maxM, ddP.a_burst, 
			 S_STABLE|S_UVTABLE|S_FLOAT| S_QUITONBOUND);
    S_tag(ddC.a_burst,"a_burst, burstiness PYP");
    if ( verbose ) S_report(ddC.a_burst,NULL);
  }
}

void cache_free() {
  if ( ddP.a_mu!=0 ) S_free(ddC.a_mu);
  if ( ddP.a_phi0!=0 ) S_free(ddC.a_phi0);
  if ( ddP.a_phi1!=0 ) S_free(ddC.a_phi1);
  if ( ddP.a_theta!=0 ) S_free(ddC.a_theta);
  if ( PCTL_BURSTY() && ddP.a_burst!=0 )
    S_free(ddC.a_burst);
  if ( ddC.a_zero )
    S_free(ddC.a_zero);
}

void cache_update(char *par) {
  if ( strcmp(par,"am")==0 && ddP.a_mu!=0 && ddC.a_mu->a!=ddP.a_mu )
    S_remake(ddC.a_mu,ddP.a_mu);
  if ( strcmp(par,"at")==0 && ddP.a_theta!=0 && ddC.a_theta->a!=ddP.a_theta )
    S_remake(ddC.a_theta,ddP.a_theta);
  if ( strcmp(par,"ap0")==0 && ddP.a_phi0!=0 && ddC.a_phi0->a!=ddP.a_phi0 )
    S_remake(ddC.a_phi0,ddP.a_phi0);
  if ( strcmp(par,"ap1")==0 && ddP.a_phi1!=0 && ddC.a_phi1->a!=ddP.a_phi1 )
    S_remake(ddC.a_phi1,ddP.a_phi1);
  if ( PCTL_BURSTY() && strcmp(par,"ab")==0 && ddC.a_burst->a!=ddP.a_burst )
    S_remake(ddC.a_burst,ddP.a_burst);
}
