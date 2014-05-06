/*
 * Likelihood calculations
 * Copyright (C) 2013-2014 Jinjing Li and Wray Buntine
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Jinjing Li (jinjingli@gmail.com)
 *         Wray Buntine (wray.buntine@monash.edu)
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "yap.h"
#include "util.h"
#include "stable.h"
#include "lgamma.h"
#include "tca.h"
#include "stats.h"

/*
 *     This is the central likelihood function;
 *
 *     * effects appear in likesub.c
 *     * extracts appear in sampleX.c
 *
 *     uses:   X tables, Y tables and caches
 *
 *     This means the perplexities computed from this will
 *     not calibrate with actual perplexities for test data.
 */

/*
 *     the binomial co-efficients are just a device introduced to
 *     allow sampling, so we don't usually include them in the likelihood
 */

// Pochhammer in log
double poch(double b, double a, int N) {
    int n;
    double result;
    if ( N<=0 )
      return 0.0;
    if ( a==0 ) {
      return N*log(b);
    }
    if ( N<=3 ) {
      result = 0;
      for (n=0; n<N; n++) {
        result += log(b+n*a);
      }
    } else {
      result = N*log(a) + lgamma(b/a+N) - lgamma(b/a);
    }
    return result;
}

double likelihood_PYalpha() {
  int e,d,t;
  double likelihood = 0;
  double lgb = lgamma(ddP.b_theta);

  // theta part
  likelihood += pctl_gammaprior(ddP.b_theta);
  for (d=0; d<ddN.DT; d++) {
    likelihood += poch(ddP.b_theta, ddP.a_theta, ddS.C_dT[d]);
    likelihood -= gammadiff(ddS.N_dT[d], ddP.b_theta, lgb);
    for (t=0; t<ddN.T; t++) {
      if ( ddS.n_dt[d][t]>1 ) 
	likelihood += S_S(ddC.a_theta, ddS.n_dt[d][t], ddS.c_dt[d][t]);
    }
    yap_infinite(likelihood);
  }

  // mu0 part, use a Dirichlet
  for (t=0; t<ddN.T; t++) {
    likelihood += gammadiff(ddS.cp_et[0][t], ddP.b_mu0/ddN.T, 0);
  }
  likelihood -= gammadiff(ddS.Cp_e[0], ddP.b_mu0, 0);
  yap_infinite(likelihood);

  // mu part
  likelihood += pctl_gammaprior(ddP.b_mu[0]);
  if ( ddN.E>1 )
    likelihood += pctl_gammaprior(ddP.b_mu[1]);
  for (e=0; e<ddN.E; e++) {
    likelihood += poch(ddP.b_mu[e], ddP.a_mu, ddS.Cp_e[e]);
    if (e<ddN.E-1) {
      likelihood -= gammadiff(ddS.C_e[e]+ddS.Cp_e[e+1], ddP.b_mu[e], 0);
    } else {
      likelihood -= gammadiff(ddS.C_e[e], ddP.b_mu[e], 0);
    }
    for (t=0; t<ddN.T; t++) {
      if ( ddS.cp_et[e][t]==0 )
	continue;
      if (e==ddN.E-1) {
      	// last epoch
      	likelihood += S_S(ddC.a_mu, ddS.C_eDt[e][t], ddS.cp_et[e][t]);
      } else {
      	likelihood += S_S(ddC.a_mu, ddS.C_eDt[e][t] + ddS.cp_et[e+1][t], ddS.cp_et[e][t]);
      }
    }
    yap_infinite(likelihood);
  }
  
  return likelihood;
}

double likelihood_PYbeta() {
  int e,t,v;
  double likelihood = 0;

  // prior
  for (t=0; t<ddN.T; t++) {
    likelihood += pctl_gammaprior(ddP.b_phi[0][t]);
    if ( ddN.E>1 )
      likelihood += pctl_gammaprior(ddP.b_phi[1][t]);
  }
  //  note we dont sample a_phi and b_phi0
  
  // phi0 part
  likelihood += poch(ddP.b_phi0, ddP.a_phi, ddS.S_0_nz);
  likelihood -= gammadiff(ddS.S_0, ddP.b_phi0, 0);
  for (v=0; v<ddN.W; v++) {
    if ( ddS.S_0vT[v]>0 )
      likelihood += S_S(ddC.a_phi, ddS.S_0vT[v], 1);
  }
  yap_infinite(likelihood);

  // phi part
  for (e=0; e<ddN.E; e++) {
    for (t=0; t<ddN.T; t++) {
      if ( ddS.S_eVt[e][t]==0 )
	continue;
      likelihood += poch(ddP.b_phi[e][t], ddP.a_phi, ddS.S_eVt[e][t]);
      if (e<ddN.E-1) {
	likelihood -=
	  gammadiff(ddS.M_eVt[e][t] + ddS.S_eVt[e+1][t], ddP.b_phi[e][t], 0);
      } else {
	likelihood -= gammadiff(ddS.M_eVt[e][t], ddP.b_phi[e][t], 0);
      }
      yap_infinite(likelihood);
      for (v=0; v<ddN.W; v++) {
	if ( ddS.s_evt[e][v][t]==0 )
	  continue;
	if (e<ddN.E-1) {
	  likelihood += S_S(ddC.a_phi, ddS.m_evt[e][v][t] + ddS.s_evt[e+1][v][t] , ddS.s_evt[e][v][t]);
	} else {
	  likelihood += S_S(ddC.a_phi, ddS.m_evt[e][v][t], ddS.s_evt[e][v][t]);
	}
      }
      yap_infinite(likelihood);
    }
  }

  return likelihood;
}

double likelihood_burst() {
  double b[ddN.T];
  int t;
  for (t=0; t<ddN.T; t++)
    b[t] = ddP.b_burst;
  return dmi_likelihood(&ddM,pctl_gammaprior,ddP.a_burst,b,ddC.a_burst);
}

double likelihood() {
  double likelihood = 0;
  /*
   *  doc X topic part
   */
  likelihood += likelihood_PYalpha();

  /*
   *  word X topic part
   */
  likelihood += likelihood_PYbeta();

  /*
   *  doc burstiness
   */
  if ( PCTL_BURSTY() )
    likelihood += likelihood_burst();

  yap_infinite(likelihood);
  return likelihood;
}


