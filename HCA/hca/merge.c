/*
 * Likelihood diff due to merge
 * Copyright (C) 2014 Wray Buntine
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
 *    adapted by diffing code in like.c
 */ 

double likemerge_bdk(k1,k2) {
  yap_quit("NI\n");
  // return dmi_likemerge(&ddM,pctl_gammaprior,ddP.ad, ddP.bdk,ddC.SD);
}

double likemerge_DIRalpha(k1,k2) {
  /*
   *   Dirichlet for topics
   */
  int i,t;
  double likelihood = 0;
  for (i=0; i<ddN.DT; i++) {
    if ( ddS.Ndt[i][k2]>0 ) {
      likelihood +=
	gammadiff((int)ddS.Ndt[i][k1]+ddS.Ndt[i][k2], ddP.alphapr[k1], 0.0)
	- gammadiff((int)ddS.Ndt[i][k1], ddP.alphapr[k1], 0.0)
	- gammadiff((int)ddS.Ndt[i][k2], ddP.alphapr[k2], 0.0);
    }
  }
  yap_infinite(likelihood);
  return likelihood;
}

double likemerge_PYalpha(k1,k2) {
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

double likemerge_PYalpha_HDP(k1,k2) {
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

double likemerge_PYalpha_PDP(k1,k2) {
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

double likemerge_PYalpha_HPDD(k1,k2) {
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

double likemerge_DIRbeta(k1,k2) {
  int j,t;
  double likelihood = 0;
  double val = 0;
  for (j=0; j<ddN.W; j++) {
    if ( ddS.Nwt[j][k2]>0 ) {
      assert(ddP.betapr[j]>0);
      val += gammadiff((int)ddS.Nwt[j][k1]+ddS.Nwt[j][k2], ddP.betapr[j], 0.0)
	- gammadiff((int)ddS.Nwt[j][k1], ddP.betapr[j], 0.0)
	- gammadiff((int)ddS.Nwt[j][k2], ddP.betapr[j], 0.0);
    } 
  }     
  val -= gammadiff((int)ddS.NWt[k1]+ddS.NWt[k2], ddP.betatot, 0.0)
    - gammadiff((int)ddS.NWt[k1], ddP.betatot, 0.0)
    - gammadiff((int)ddS.NWt[k2], ddP.betatot, 0.0);
  likelihood += val;
  yap_infinite(likelihood);
  return likelihood;
}

double likemerge_PYbeta(k1,k2) {
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

double likemerge_PYbeta_PDP(k1,k2) {
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

double likemerge_PYbeta_HDP(k1,k2) {
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

double likemerge_PYbeta_HPDD(k1,k2) {
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

/*
 *  place all topic k2 data in topic k1
 */
double likemerge(int k1, int k2) {
  int t,j;
  double likelihood = 0;
  if ( k1<=0 || k2<=0 ) 
    return 0.0;
  /*
   *   PYP doc part
   */
  if ( ddP.bdk!=NULL ) 
    likelihood += likemerge_bdk(k1,k2);

  /*
   *  doc X topic part
   */
  if ( ddP.PYalpha ) {
    likelihood += likemerge_PYalpha(k1,k2);
    /*
     *  term for root node
     */
    if ( ddP.PYalpha==H_PDP )
      likelihood += likemerge_PYalpha_PDP(k1,k2);
    else if ( ddP.PYalpha==H_HDP ) 
      likelihood += likemerge_PYalpha_HDP(k1,k2);
    else 
      likelihood += likemerge_PYalpha_HPDD(k1,k2);
  } else 
    likelihood += likemerge_DIRalpha(k1,k2);
  
  /*
   *  word X topic part
   */
  if ( ddP.phi!=NULL ) {
    /*
     *    no learning, just preexisting phi[][]
     */
    for (j=0; j<ddN.W; j++) 
      if ( ddS.Nwt[j][k2]>0 ) 
	likelihood += ddS.Nwt[j][k2]*(log(ddP.phi[k1][j])-log(ddP.phi[k2][j]));
  } else if ( ddP.PYbeta ) {
    likelihood += likemerge_PYbeta(k1,k2);
    /*
     *  term for root node
     */
    if ( ddP.PYbeta==H_PDP ) 
      likelihood += likemerge_PYbeta_PDP(k1,k2);
    else if ( ddP.PYbeta==H_HDP ) 
      likelihood += likemerge_PYbeta_HDP(k1,k2);
    else 
      likelihood += likemerge_PYbeta_HPDD(k1,k2);
  } else 
    likelihood += likemerge_DIRbeta(k1,k2);
 
  yap_infinite(likelihood);
  return likelihood;
}


