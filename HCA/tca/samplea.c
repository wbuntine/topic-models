/*
 * Sampling utility for a
 * Copyright (C) 2013 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@monash.edu)
 *
 */


#include <stdio.h>  
#include <unistd.h>   
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "yap.h"
#include "lgamma.h"
#include "myarms.h"
#include "stable.h"
#include "sample.h"
#include "stats.h"

/*
 *  debugging uses static memory, so only allow without threads
 */
#ifndef H_THREADS
// #define A_DEBUG
#endif

#ifdef A_DEBUG
  static double last_val = 0;
  static double last_like = 0;
#endif
/*
 *  globals from tca.h and stats.h
 */
extern int verbose;

/*
 *  defined in like.c
 */
double poch(double b, double a, int N);

/*
 */
static double aterms_theta(double mya, void *mydata) {
  int i, t;
  double val = 0;
#ifdef A_DEBUG
  float save_a = ddC.a_theta->a;
  double like;
#endif
  S_remake(ddC.a_theta, mya);
  for (i=0; i<ddN.DT; i++) {
    for (t=0; t<ddN.T; t++) {
      if ( ddS.n_dt[i][t]>1 ) {
	val += S_S(ddC.a_theta,ddS.n_dt[i][t],ddS.c_dt[i][t]);
      }
    }
    val += poch(ddP.b_theta, mya, ddS.C_dT[i]);
  }  
  myarms_evals++;
#ifdef A_DEBUG
  yap_message("Eval aterms_theta(%lf) = %lf (S had %f)", mya, val, save_a);
  ddP.a_theta = mya;
  cache_update("at");
  like = likelihood();
  if ( last_val != 0 ) {
    yap_message(", lp=%lf diffs=%lf vs %lf\n", like, 
		val-last_val, like-last_like);
  }
  last_like = like;
  last_val = val;
#endif
  return val;
}

static double aterms_mu(double mya, void *mydata) {
  int e, t;
  double val = 0;
#ifdef A_DEBUG
  float save_a = ddC.a_mu->a;
  double like;
#endif
  S_remake(ddC.a_mu, mya);
  for (e=0; e<ddN.E; e++) {
    for (t=0; t<ddN.T; t++) {
      if ( ddS.cp_et[e][t]==0 )
	continue;
      if (e==ddN.E-1) 
	val += S_S(ddC.a_mu, ddS.C_eDt[e][t], ddS.cp_et[e][t]); 
      else
	val += 
	  S_S(ddC.a_mu, ddS.C_eDt[e][t] + ddS.cp_et[e+1][t], ddS.cp_et[e][t]);
    }
    val += poch(ddP.b_mu[e], mya, ddS.Cp_e[e]);
  }  
  myarms_evals++;
#ifdef A_DEBUG
  yap_message("Eval aterms_mu(%lf) = %lf (S had %f)", mya, val, save_a);
  ddP.a_mu = mya;
  cache_update("am");
  like = likelihood();
  if ( last_val != 0 ) {
    yap_message(", lp=%lf diffs=%lf vs %lf\n", like, 
		val-last_val, like-last_like);
  }
  last_like = like;
  last_val = val;
#endif
  return val;
}

/*
 */
static double aterms_burst(double mya, void *mydata) {
  double b[ddM.T];
  double val = 0;
  uint16_t **docstats = (uint16_t **)mydata;
#ifdef A_DEBUG
  float save_a = ddC.a_burst->a;
  double like;
#endif
  int t;
  for (t=0; t<ddM.T; t++)
    b[t] = ddP.b_burst;
  cache_update("ab");
  val = dmi_likelihood_aterms(&ddM, docstats,
			      pctl_gammaprior, mya, b, ddC.a_burst);
  myarms_evals++;
#ifdef A_DEBUG
  yap_message("Eval adkterms(%lf) = %lf", mya, val);
  like = dmi_likelihood(&ddM,pctl_gammaprior,mya,b,ddC.a_burst);
  if ( last_val != 0 ) {
    yap_message(", lp=%lf diffs=%lf vs %lf\n", like, 
		val-last_val, like-last_like);
  } else
    yap_message("\n");
  last_like = like;
  last_val = val;
#endif
  return val;
}


static double aterms_phi0(double mya, void *mydata) {
  int v;
  double val = 0;
#ifdef A_DEBUG
  float save_a = ddC.a_phi0->a;
  double like;
#endif
  S_remake(ddC.a_phi0, mya);
  val += poch(ddP.b_phi0, mya, ddS.S_0_nz);
  for (v=0; v<ddN.W; v++) {
    if ( ddS.S_0vT[v]>1 )
      val += S_S(ddC.a_phi0, ddS.S_0vT[v], 1);
  }
  myarms_evals++;
#ifdef A_DEBUG
  yap_message("Eval aterms_phi0(%lf) = %lf (S had %f)", mya, val, save_a);
  ddP.a_phi0 = mya;
  cache_update("ap0");
  like = likelihood();
  if ( last_val != 0 ) {
    yap_message(", lp=%lf diffs=%lf vs %lf\n", like, 
		val-last_val, like-last_like);
  }
  last_like = like;
  last_val = val;
#endif
  return val;
}
static double aterms_phi1(double mya, void *mydata) {
  int e, t, v;
  double val = 0;
#ifdef A_DEBUG
  float save_a = ddC.a_phi1->a;
  double like;
#endif
  S_remake(ddC.a_phi1, mya);
  for (e=0; e<ddN.E; e++) {
    for (t=0; t<ddN.T; t++) {
      if ( ddS.S_Vte[t][e]==0 )
	continue;
      val += poch(ddP.b_phi[e][t], mya, ddS.S_Vte[t][e]);
      for (v=0; v<ddN.W; v++) {
	if ( ddS.s_vte[v][t][e]==0 )
	  continue;
	if (e<ddN.E-1) {
	  val += S_S(ddC.a_phi1, ddS.m_vte[v][t][e] + ddS.s_vte[v][t][e+1] , ddS.s_vte[v][t][e]);
	} else {
	  val += S_S(ddC.a_phi1, ddS.m_vte[v][t][e], ddS.s_vte[v][t][e]);
	}
      }
    }
  }
  myarms_evals++;
#ifdef A_DEBUG
  yap_message("Eval aterms_phi1(%lf) = %lf (S had %f)", mya, val, save_a);
  ddP.a_phi1 = mya;
  cache_update("ap"1);
  like = likelihood();
  if ( last_val != 0 ) {
    yap_message(", lp=%lf diffs=%lf vs %lf\n", like, 
		val-last_val, like-last_like);
  }
  last_like = like;
  last_val = val;
#endif
  return val;
}


void sample_at(double *mya) {
#ifdef A_DEBUG
  last_val = 0;
  last_like = 0;
#endif
  myarms(PYP_DISC_MIN, PYP_DISC_MAX, &aterms_theta, NULL, mya, "at");
  ddP.a_theta = *mya;
  cache_update("at");
}

void sample_am(double *mya) {
#ifdef A_DEBUG
  last_val = 0;
  last_like = 0;
#endif
  myarms(PYP_DISC_MIN, PYP_DISC_MAX, &aterms_mu, NULL, mya, "am");
  ddP.a_mu = *mya;
  cache_update("am");
 }

void sample_ab(double *mya) {
  uint16_t **docstats;
#ifdef A_DEBUG
  last_val = 0;
  last_like = 0;
#endif  
  docstats = dmi_astore(&ddM);
  myarms(PYP_DISC_MIN, PYP_DISC_MAX, &aterms_burst, docstats, mya, "ab");
  ddP.a_burst = *mya;
  cache_update("ab");
  dmi_freeastore(&ddM, docstats);
}

void sample_ap1(double *mya) {
#ifdef A_DEBUG
  last_val = 0;
  last_like = 0;
#endif
  myarms(PYP_DISC_MIN, PYP_DISC_MAX, &aterms_phi1, NULL, mya, "ap1");
  ddP.a_phi1 = *mya;
  cache_update("ap1");
}
void sample_ap0(double *mya) {
#ifdef A_DEBUG
  last_val = 0;
  last_like = 0;
#endif
  myarms(PYP_DISC_MIN, PYP_DISC_MAX, &aterms_phi0, NULL, mya, "ap0");
  ddP.a_phi0 = *mya;
  cache_update("ap0");
}


