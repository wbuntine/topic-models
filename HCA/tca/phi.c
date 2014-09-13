/*
 * Computing, updating and saving topicXword/phi estimates
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
 *  If ddP.memory is set, then statistics kept
 *  in file in binary format.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "yap.h"
#include "util.h"
#include "tca.h"
#include "data.h"
#include "stats.h"

static char *phi_file = NULL;
static char *mu_file = NULL;

void phi_init(char *resstem) {
  int e;
  phi_file = yap_makename(resstem,".phi");
  ddS.phi_cnt = 0;
  ddS.phi = malloc(sizeof(ddS.phi[0])*ddN.E);
  for (e=0; e<ddN.E; e++) {
    ddS.phi[e] = fmat(ddN.W, ddN.T);
    if ( !ddS.phi[e] )
      yap_quit("Not enough memory in phi_init()\n");
  }
}
void mu_init(char *resstem) {
  mu_file = yap_makename(resstem,".mu");
  ddS.mu_cnt = 0;
  ddS.mu = fmat(ddN.E,ddN.T);
  if ( !ddS.mu )
    yap_quit("Not enough memory in mu_init()\n");
}

void phi_free() {
  ddS.phi_cnt = 0;
  if ( phi_file ) {
    free(phi_file);
    phi_file = NULL;
  }
  if ( ddS.phi ) {
    int e;  
    for (e=0; e<ddN.E; e++) {
      free(ddS.phi[e][0]); free(ddS.phi[e]);
    }
    free(ddS.phi);
    ddS.phi = NULL;
  }
}

void mu_free() {
  ddS.mu_cnt = 0;
  if ( mu_file ) {
    free(mu_file);
    mu_file = NULL;
  }
  if ( ddS.mu ) {
    free(ddS.mu[0]);free(ddS.mu[0]);
    ddS.mu = NULL;
  }
}

/*
 *    only need to save if *not* storing in file
 */
void phi_save() {
  int e;
  char ebuf[10];
  char *fname;
  for (e=0; e<ddN.E; e++) {
    sprintf(ebuf,"%03d",e);
    fname = yap_makename(phi_file, ebuf);
    write_fmat(fname,ddN.W,ddN.T,ddS.phi[e]);
    free(fname);
  }
}

void mu_save() {
  write_fmat(mu_file,ddN.E,ddN.T,ddS.mu);
}

void phi_update() {
  int e, t, w;

  ddS.phi_cnt++;
}

void mu_update() {
  int e, t;
  ;
  ddS.mu_cnt++;
}
