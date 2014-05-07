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

// #define A_DEBUG

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

static uint16_t **docstats;
/*
 */
static double aterms_burst(double mya, void *mydata) {
  double b[ddM.T];
  double val = 0;
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


static double aterms_phi(double mya, void *mydata) {
  int e, t, v;
  double val = 0;
#ifdef A_DEBUG
  float save_a = ddC.a_phi->a;
  double like;
#endif
  S_remake(ddC.a_phi, mya);
  val += poch(ddP.b_phi0, mya, ddS.S_0_nz);
  for (v=0; v<ddN.W; v++) {
    if ( ddS.S_0vT[v]>1 )
      val += S_S(ddC.a_phi, ddS.S_0vT[v], 1);
  }
  for (e=0; e<ddN.E; e++) {
    for (t=0; t<ddN.T; t++) {
      if ( ddS.S_eVt[e][t]==0 )
	continue;
      val += poch(ddP.b_phi[e][t], mya, ddS.S_eVt[e][t]);
      for (v=0; v<ddN.W; v++) {
	if ( ddS.s_evt[e][v][t]==0 )
	  continue;
	if (e<ddN.E-1) {
	  val += S_S(ddC.a_phi, ddS.m_evt[e][v][t] + ddS.s_evt[e+1][v][t] , ddS.s_evt[e][v][t]);
	} else {
	  val += S_S(ddC.a_phi, ddS.m_evt[e][v][t], ddS.s_evt[e][v][t]);
	}
      }
    }
  }
  myarms_evals++;
#ifdef A_DEBUG
  yap_message("Eval aterms_phi(%lf) = %lf (S had %f)", mya, val, save_a);
  ddP.a_phi = mya;
  cache_update("ap");
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
  if ( verbose>1 )
    yap_message("sample_at (pre):  a_theta=%lf, lp=%lf\n",
		*mya, likelihood());
  myarms(PYP_DISC_MIN, PYP_DISC_MAX, &aterms_theta, NULL, mya, "at");
  ddP.a_theta = *mya;
  cache_update("at");
  if ( verbose>1 )
    yap_message("sample_at (post):  a_theta=%lf, lp=%lf\n",
		*mya, likelihood());
}

void sample_am(double *mya) {
#ifdef A_DEBUG
  last_val = 0;
  last_like = 0;
#endif
  if ( verbose>1 )
    yap_message("sample_am (pre):  a_mu=%lf, lp=%lf\n",
		ddP.a_mu, likelihood());
  myarms(PYP_DISC_MIN, PYP_DISC_MAX, &aterms_mu, NULL, mya, "am");
  ddP.a_mu = *mya;
  cache_update("am");
  if ( verbose>1 )
    yap_message("sample_am (post):  a_mu=%lf, lp=%lf\n",
		ddP.a_mu, likelihood());
 }

void sample_ab(double *mya) {
#ifdef A_DEBUG
  last_val = 0;
  last_like = 0;
#endif  
  docstats = dmi_astore(&ddM);
  if ( verbose>1 )
    yap_message("sample_ab (pre):  a_burst=%lf, lp=%lf\n",
		*mya, likelihood());
  myarms(PYP_DISC_MIN, PYP_DISC_MAX, &aterms_burst, NULL, mya, "ab");
  ddP.a_burst = *mya;
  cache_update("ab");
  dmi_freeastore(&ddM, docstats);
  if ( verbose>1 )
    yap_message("sample_ab (post):  a_burst=%lf, lp=%lf\n",
		*mya, likelihood());
}

void sample_ap(double *mya) {
#ifdef A_DEBUG
  last_val = 0;
  last_like = 0;
#endif
  if ( verbose>1 )
    yap_message("sample_ap (pre):  a_phi=%lf, lp=%lf\n",
		*mya, likelihood());
  myarms(PYP_DISC_MIN, PYP_DISC_MAX, &aterms_phi, NULL, mya, "ap");
  ddP.a_phi = *mya;
  cache_update("ap");
  if ( verbose>1 )
    yap_message("sample_ap (post):  a_phi=%lf, lp=%lf\n",
		*mya, likelihood());
}


