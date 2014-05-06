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

void cache_init() {
  ddC.a_mu = S_make(ddP.maxN/10, 100, ddP.maxN+1, ddP.maxM, ddP.a_mu, 
#ifdef H_THREADS
 		    S_THREADS|
#endif
		    S_STABLE|S_UVTABLE|S_FLOAT| S_QUITONBOUND);
  ddC.a_theta = S_make(ddD.N_dTmax/2, 100, ddD.N_dTmax+1, ddP.maxM, ddP.a_theta,
#ifdef H_THREADS
                       S_THREADS|
#endif
		       S_STABLE|S_UVTABLE|S_FLOAT| S_QUITONBOUND);
  ddC.a_phi = S_make(ddP.maxN/10, ddP.maxM/2, ddP.maxN+1, ddP.maxM, ddP.a_phi, 
#ifdef H_THREADS
                     S_THREADS|
#endif
                     S_STABLE|S_UVTABLE|S_FLOAT| S_QUITONBOUND);
  if ( !PCTL_BURSTY() )
    ddC.a_burst = NULL;
  else
    ddC.a_burst = S_make(ddD.N_dTmax/2, 100, ddD.N_dTmax+1, ddP.maxM, ddP.a_burst, 
			 S_STABLE|S_UVTABLE|S_FLOAT| S_QUITONBOUND);
  S_tag(ddC.a_mu,"a_mu, topic time PYP");
  S_tag(ddC.a_theta,"a_theta, topic->doc PYP");
  S_tag(ddC.a_phi,"a_phi, word time PYP");
  if ( PCTL_BURSTY() )
    S_tag(ddC.a_burst,"a_burst, burstiness PYP");
  /*
   *   now the libstb library has its own separate error handling,
   *   so make sure it goes through our yapper
   */
  yaps_yapper(yap_va);
  if ( verbose ) {
    S_report(ddC.a_mu,NULL);
    S_report(ddC.a_theta,NULL);
    S_report(ddC.a_phi,NULL);
    if ( PCTL_BURSTY() )
      S_report(ddC.a_burst,NULL);
  }
}

void cache_free() {
  S_free(ddC.a_mu);
  S_free(ddC.a_phi);
  S_free(ddC.a_theta);
  if ( PCTL_BURSTY() )
    S_free(ddC.a_burst);
}

void cache_update(char *par) {
  if ( strcmp(par,"am")==0 && ddC.a_mu->a!=ddP.a_mu )
    S_remake(ddC.a_mu,ddP.a_mu);
  if ( strcmp(par,"at")==0 && ddC.a_theta->a!=ddP.a_theta )
    S_remake(ddC.a_theta,ddP.a_theta);
  if ( strcmp(par,"ap")==0 && ddC.a_phi->a!=ddP.a_phi )
    S_remake(ddC.a_phi,ddP.a_phi);
  if ( PCTL_BURSTY() && strcmp(par,"ab")==0 && ddC.a_burst->a!=ddP.a_burst )
    S_remake(ddC.a_burst,ddP.a_burst);
}
