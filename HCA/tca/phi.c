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
static char *theta_file = NULL;
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
void theta_init(char *resstem) {
  theta_file = yap_makename(resstem,".theta");
  ddS.theta = fmat(ddN.D, ddN.T);
  if ( !ddS.theta )
    yap_quit("Not enough memory in theta_init()\n");
  ddS.theta_cnt = 0;
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

void theta_free() {
  ddS.theta_cnt = 0;
  if ( theta_file ) {
    free(theta_file);
    theta_file = NULL;
  }
  if ( ddS.theta ) {
    free(ddS.theta[0]); free(ddS.theta);
    ddS.theta = NULL;
  }
}

void mu_free() {
  ddS.mu_cnt = 0;
  if ( mu_file ) {
    free(mu_file);
    mu_file = NULL;
  }
  if ( ddS.mu ) {
    free(ddS.mu[0]);free(ddS.mu);
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
  if ( verbose ) 
    yap_message("Saving '%s' with count %d\n", phi_file, ddS.phi_cnt);
  for (e=0; e<ddN.E; e++) {
    sprintf(ebuf,"%03d",e);
    fname = yap_makename(phi_file, ebuf);
    write_fmat(fname,ddN.W,ddN.T,ddS.phi[e]);
    free(fname);
  }
}
void theta_save() {
  if ( verbose ) 
    yap_message("Saving '%s' with count %d\n", theta_file, ddS.theta_cnt);
  write_fmat(theta_file,ddN.D,ddN.T,ddS.theta);
}

void mu_save() {
  if ( verbose ) 
    yap_message("Saving '%s' with count %d\n", mu_file, ddS.mu_cnt);
  write_fmat(mu_file,ddN.E,ddN.T,ddS.mu);
}

void phi_update() {
  int e, t, w;
  double **wmtx;
  
  wmtx = dmat(ddN.W, ddN.T);
  phi_prob_iter(-1, wmtx);
  for (e=0; e<ddN.E; e++) {
    phi_prob_iter(e, wmtx);
    for (w=0; w<ddN.W; w++)
      for (t=0; t<ddN.T; t++)
        ddS.phi[e][w][t] = (ddS.phi_cnt*ddS.phi[e][w][t] + wmtx[w][t]) 
          / (ddS.phi_cnt+1);
  }
  free(wmtx[0]); free(wmtx);
  ddS.phi_cnt++;
}

void mu_update() {
  int e, t;
  double *pvec = dvec(ddN.T);
  mu_prob_iter(-1, pvec);

  for (e=0; e<ddN.E; e++) {
    mu_prob_iter(e, pvec);
    for (t=0; t<ddN.T; t++)
      ddS.mu[e][t] = (ddS.mu_cnt*ddS.mu[e][t] + pvec[t]) / (ddS.mu_cnt+1);
  }
  free(pvec);
  ddS.mu_cnt++;
}

void theta_update() {
  int d, e, t;
  double *pvec = dvec(ddN.T);
  mu_prob_iter(-1, pvec);

  for (d=0,e=0; e<ddN.E; e++) {
    mu_prob_iter(e, pvec);
    for ( ; d<ddN.DT && ddD.e[d]==e; d++) {
      double total = ddS.N_dT[d] + ddP.b_theta;
	for (t=0; t<ddN.T; t++) {
	  double prob = (ddS.n_dt[d][t]-ddP.a_theta*ddS.c_dt[d][t]) 
	    + (ddP.b_theta + ddP.a_theta * ddS.C_dT[d])*pvec[t];
	  prob /= total;
	  ddS.theta[d][t] = (ddS.theta_cnt*ddS.theta[d][t] + prob) 
	    / (ddS.theta_cnt+1);
	}
    }
  }
  if ( ddN.D>ddN.DT ) {
    for (e=0; e<ddN.E; e++) {
      mu_prob_iter(e, pvec);
      for ( ; d<ddN.D && ddD.e[d]==e; d++) {
	double total = ddS.N_dT[d] + ddP.b_theta;
	for (t=0; t<ddN.T; t++) {
	  double prob = (ddS.n_dt[d][t]-ddP.a_theta*ddS.c_dt[d][t]) 
	    + (ddP.b_theta + ddP.a_theta * ddS.C_dT[d])*pvec[t];
	  prob /= total;
	  ddS.theta[d][t] = (ddS.theta_cnt*ddS.theta[d][t] + prob) 
	    / (ddS.theta_cnt+1);
	}
      }
    }
  }
  free(pvec);
  ddS.theta_cnt++;
}

float *phi_mean(int k) {
  int e, w;
  float tot=0, prop[ddN.E];
  float *wvec;
  if ( !ddS.phi )
    return NULL;
  wvec = fvec(ddN.W);
  for (e=0; e<ddN.E; e++)
    tot += prop[e] = ddD.esize[e];
  for (e=0; e<ddN.E; e++)
    prop[e] /= tot;
  for (e=0; e<ddN.E; e++)
    for (w=0; w<ddN.W; w++)
      wvec[w] += ddS.phi[e][w][k] * prop[e];
  return wvec;
}

float *mu_mean() {
  int e, k;
  float tot=0, prop[ddN.E];
  float *kvec;
  if ( !ddS.mu )
    return NULL;
  kvec = fvec(ddN.T);
  for (e=0; e<ddN.E; e++)
    tot += prop[e] = ddD.esize[e];
  for (e=0; e<ddN.E; e++)
    prop[e] /= tot;
  for (e=0; e<ddN.E; e++)
    for (k=0; k<ddN.T; k++)
      kvec[k] += ddS.mu[e][k] * prop[e];
  return kvec;
}

float *theta_mean() {
  int d, k;
  float *kvec;
  if ( !ddS.theta )
    return NULL;
  kvec = fvec(ddN.T);
  for (d=0; d<ddN.D; d++)
    for (k=0; k<ddN.T; k++)
      kvec[k] += ddS.theta[d][k];
  for (k=0; k<ddN.T; k++)
    kvec[k] /= ddN.D;
  return kvec;
}
