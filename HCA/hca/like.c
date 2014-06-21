/*
 * Likelihood calculations
 * Copyright (C) 2009-2014 Wray Buntine
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@monash.edu)
 *
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
#include "hca.h"
#include "stats.h"
#include "cache.h"

/*
 *     this is the central likelihood function;
 *
 *     * effects appear in likesub.c
 *     * extracts appear in sampleX.c
 *
 *     uses:   X tables, Y tables and caches
 *             lgba, lgb, lgalphac, lgalphatot, lgbaw, lgbw, lgbetac
 *             lgbetatot
 *
 *     NB.    No prior terms for hyperparameters used.
 *
 *            No prior used for alpha, beta, gamma0/1.
 *            Uniform prior valid for a, a0, aw, aw0
 *            bw, bw0 have an inverse (scale) prior in samplers,
 *              but its not included in likelihood.
 *           
 *     This means the perplexities computed from this will
 *     not calibrate with actual perplexities for test data.
 */ 

#define L_CACHE

double likelihood_bdk() {
  return dmi_likelihood(&ddM,pctl_gammaprior,ddP.ad, ddP.bdk,ddC.SD);
}

double likelihood_DIRalpha() {
  /*
   *   Dirichlet for topics
   */
  int i,t;
  double likelihood = 0;
  for (i=0; i<ddN.DT; i++) {
    for (t=0; t<ddN.T; t++) {
      if ( ddS.Ndt[i][t]>0 ) {
#ifdef L_CACHE
        if ( ddP.alphac>0 )
          likelihood += gcache_value(&ddC.lgalphac, (int)ddS.Ndt[i][t]);
        else
#endif
          likelihood += gammadiff((int)ddS.Ndt[i][t], ddP.alphapr[t], 0.0);
      }
    }
#ifdef L_CACHE
    likelihood -= gcache_value(&ddC.lgalphatot, (int)ddS.NdT[i]);
#else
    likelihood -= gammadiff((int)ddS.NdT[i], ddP.alphatot, 0.0);
#endif
  }
  yap_infinite(likelihood);
  return likelihood;
}

double likelihood_PYalpha() {
  int i,t;
  double likelihood = 0;
  double la = 0;
  double lb = log(ddP.bpar);
  if ( ddP.apar>0 ) la = log(ddP.apar);
  for (i=0; i<ddN.DT; i++) {
    uint16_t Td_ = 0;
    for (t=0; t<ddN.T; t++) {
      if ( ddS.Ndt[i][t]>0 ) {
	Td_ += ddS.Tdt[i][t];
	if ( ddS.Ndt[i][t]>1 ) {
	  likelihood += S_S(ddC.SX,ddS.Ndt[i][t],ddS.Tdt[i][t]);
	}
      }
    }
    yap_infinite(likelihood);
    if ( ddP.apar==0 ) {
      likelihood += Td_*lb;
    } else {
      likelihood += Td_*la + gcache_value(&ddC.lgba, (int)Td_);
    }
    likelihood -= gcache_value(&ddC.lgb, (int)ddS.NdT[i]);
    yap_infinite(likelihood);
  }  
  return likelihood;
}
double likelihood_PYalpha_HDP() {
  /*
   *    the DP prior, a0==0, its a Dirichlet
   */
  int t;
  double likelihood = 0;
  assert(ddP.a0==0);
  for (t=0; t<ddN.T; t++) {
    if ( ddS.TDt[t]>0 ) {
      likelihood += gammadiff((int)ddS.TDt[t], ddP.b0*ddP.alphapr[t], 0.0);
    }
  }      
  likelihood -= lgamma(ddP.b0+ddS.TDT) - lgamma(ddP.b0);
  yap_infinite(likelihood);
  return likelihood;
}

double likelihood_PYalpha_PDP() {
  /*
   *    the constant prior
   */
  int t;
  double likelihood = 0;
  for (t=0; t<ddN.T; t++) {
    if ( ddS.TDt[t]>0 ) {
      assert(ddP.alphapr[t]>0);
      likelihood += ddS.TDt[t]*log(ddP.alphapr[t]);
    }
  }      
  yap_infinite(likelihood);
  return likelihood;
}

double likelihood_PYalpha_HPDD() {
  /*
   *    the PDD prior
   */
  int t;
  double likelihood = 0;
  double l1a0 = log(1-ddP.a0);
  double l2a0 = log((1-ddP.a0)*(2-ddP.a0));
  double lga0 = lgamma(1-ddP.a0);
  /*
   *    the PDD prior
   */
  if ( ddP.a0==0 )
    likelihood += ddS.TDTnz*log(ddP.b0);
  else
    likelihood += ddS.TDTnz*log(ddP.a0) + lgamma(ddP.b0/ddP.a0+ddS.TDTnz)
      - lgamma(ddP.b0/ddP.a0);
  yap_infinite(likelihood);
  for (t=0; t<ddN.T; t++) { 
    /*  note the root node is a PDD so all t's = 1 */
    if ( ddS.TDt[t]>1 ) {
      if ( ddS.TDt[t]==2 )
	likelihood += l1a0;
      else if ( ddS.TDt[t]==3 )
	likelihood += l2a0;
      else
	likelihood += lgamma(ddS.TDt[t]-ddP.a0) - lga0;
    }
  }
  likelihood -= lgamma(ddP.b0+ddS.TDT) - lgamma(ddP.b0);
  yap_infinite(likelihood);
  return likelihood;
}

double likelihood_DIRbeta() {
  int j,t;
  double likelihood = 0;
  double val = 0;
  for (t=0; t<ddN.T; t++) {
    for (j=0; j<ddN.W; j++) {
      if ( ddS.Nwt[j][t]>0 ) {
	assert(ddP.betapr[j]>0);
#ifdef L_CACHE
	if ( ddP.betac>0 )
	  val += gcache_value(&ddC.lgbetac, (int)ddS.Nwt[j][t]);
	else
#endif
          val += gammadiff((int)ddS.Nwt[j][t], ddP.betapr[j], 0.0);
      }      
    }
#ifdef L_CACHE
    val -= gcache_value(&ddC.lgbetatot, (int)ddS.NWt[t]);
#else
    val -= gammadiff((int)ddS.NWt[t], ddP.betatot, 0.0);
#endif
  }
  likelihood += val;
  yap_infinite(likelihood);
  return likelihood;
}

double likelihood_PYbeta() {
  int i,t;
  double likelihood = 0;
  double lbw = log(ddP.bwpar);
  double law = log(ddP.awpar);
  likelihood += pctl_gammaprior(ddP.bwpar);
  /*
   *    term for k-th node
   */
  for (t=0; t<ddN.T; t++) {
    uint32_t Tw_ = 0;
    for (i=0; i<ddN.W; i++) {
      int tt = ddS.Twt[i][t];
      int nn = ddS.Nwt[i][t];
      if ( nn>0 ) {
        Tw_ += tt;
	likelihood += S_S(ddC.SY,nn,tt);
#if 0
	if ( !finite(likelihood) ) 
	  yap_quit("Inf:  Nwt[%d][%d]=%d  Twt[i][t]=%d S.M=%d S.N=%d\n",
		   i, t, (int)ddS.Nwt[i][t],(int)ddS.Twt[i][t],ddC.SY->usedM, ddC.SY->usedN);
#endif
      }
    }
    if ( ddP.awpar==0 ) {
      likelihood += Tw_*lbw;
    } else {
#ifdef L_CACHE
      likelihood += Tw_*law + gcache_value(&ddC.lgbaw, (int)Tw_);
#else
      likelihood += Tw_*law + gammadiff((int)Tw_, ddP.bwpar/ddP.awpar, 0.0);
#endif
    }
#ifdef L_CACHE
    likelihood -= gcache_value(&ddC.lgbw, (int)ddS.NWt[t]);
#else
    likelihood -= gammadiff((int)ddS.NWt[t], ddP.bwpar, 0.0);
#endif
  }  
  yap_infinite(likelihood);   
  return likelihood;
}

double likelihood_PYbeta_PDP() {
  /*
   *    the constant prior
   */
  int j;
  double likelihood = 0;
  for (j=0; j<ddN.W; j++) {
    if ( ddS.TwT[j]>0 ) {
      likelihood += ddS.TwT[j]*log(ddP.betapr[j]);
    }
  }      
  yap_infinite(likelihood);
  return likelihood;
}

double likelihood_PYbeta_HDP() {
  /*
   *    the DP prior, aw0==0, its a Dirichlet
   */
  int j;
  double likelihood = 0;
  assert(ddP.aw0==0);
  for (j=0; j<ddN.W; j++) {
    if ( ddS.TwT[j]>0 ) {
      double p0 = ddP.bw0*ddP.betapr[j];
      assert(p0>0);
      likelihood += gammadiff((int)ddS.TwT[j],p0,0.0);
    }
  }      
  likelihood += pctl_gammaprior(ddP.bw0);
  likelihood -= lgamma(ddP.bw0+ddS.TWT) - lgamma(ddP.bw0);
  yap_infinite(likelihood);
  return likelihood;
}

double likelihood_PYbeta_HPDD() {
  /*
   *    the PDD prior
   */
  int j;
  double likelihood = 0;
  double l1aw0 = log(1-ddP.aw0);
  double l2aw0 = log((1-ddP.aw0)*(2-ddP.aw0));
  double lgaw0 = lgamma(1-ddP.aw0);
  if ( ddP.aw0==0 )
    likelihood += ddS.TWTnz*log(ddP.bw0);
  else
    likelihood += ddS.TWTnz*log(ddP.aw0) + lgamma(ddP.bw0/ddP.aw0+ddS.TWTnz)
      - lgamma(ddP.bw0/ddP.aw0);
  yap_infinite(likelihood);      
  for (j=0; j<ddN.W; j++) { 
    /*  note the root node is a PDD so all t's = 1 */
    if ( ddS.TwT[j]>1 ) {
      if ( ddS.TwT[j]==2 )
	likelihood += l1aw0;
      else if ( ddS.TwT[j]==3 )
	likelihood += l2aw0;
      else
	likelihood += lgamma(ddS.TwT[j]-ddP.aw0) - lgaw0;
    }
  }
  likelihood += pctl_gammaprior(ddP.bw0);
  likelihood -= lgamma(ddP.bw0+ddS.TWT) - lgamma(ddP.bw0);
  yap_infinite(likelihood);
  return likelihood;
}

double likelihood() {
  int t,j;
  double likelihood = 0;
  /*
   *   PYP doc part
   */
  if ( ddP.bdk!=NULL ) 
    likelihood += likelihood_bdk();

  /*
   *  doc X topic part
   */
  if ( ddP.PYalpha ) {
    likelihood += likelihood_PYalpha();
    /*
     *  term for root node
     */
    if ( ddP.PYalpha==H_PDP )
      likelihood += likelihood_PYalpha_PDP();
    else if ( ddP.PYalpha==H_HDP ) 
      likelihood += likelihood_PYalpha_HDP();
    else 
      likelihood += likelihood_PYalpha_HPDD();
  } else 
    likelihood += likelihood_DIRalpha();
  
  /*
   *  word X topic part
   */
  if ( ddP.phi!=NULL ) {
    /*
     *    no learning, just preexisting phi[][]
     */
    for (t=0; t<ddN.T; t++) 
      for (j=0; j<ddN.W; j++) 
	if ( ddS.Nwt[j][t] )
	  likelihood += ddS.Nwt[j][t]*log(ddP.phi[t][j]);
  } else if ( ddP.PYbeta ) {
    likelihood += likelihood_PYbeta();
    /*
     *  term for root node
     */
    if ( ddP.PYbeta==H_PDP ) 
      likelihood += likelihood_PYbeta_PDP();
    else if ( ddP.PYbeta==H_HDP ) 
      likelihood += likelihood_PYbeta_HDP();
    else 
      likelihood += likelihood_PYbeta_HPDD();
  } else 
    likelihood += likelihood_DIRbeta();
 
  yap_infinite(likelihood);
  return likelihood;
}


