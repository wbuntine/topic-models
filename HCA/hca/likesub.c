/*
 * Sub-variable likelihood calculations
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
 *  Specialised probability calcs for parts of the model,
 *  nothing is changed
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
#include "check.h"
#include "cache.h"

/*************************************
 *
 *     more or less mirrored routines for the 
 *        wordXtopic versus docXtopic
 *     routines for sampling or probability estimation
 */

/*
 *     Base distribution for alpha
 *
 *     alpha_update() and get_probs() call this and
 *     expect all mass on first empty topic for H_PDD
 *        i.e., no spreading, and do not call for subsequent empty topics
 *     topic_fact() and topic_prob() will not call when empty topic for H_PDD
 */
double alphabasetopicprob(int t) {
  assert(t>=0);
  if ( ddP.PYalpha==H_PDP ) 
    return ddP.alphapr[t];
  if ( ddP.PYalpha==H_HDP ) {
    assert(ddP.a0==0);
    return  ( (double)ddS.TDt[t]+ddP.b0*ddP.alphapr[t] )
      /((double)ddS.TDT+ddP.b0);
  }
  /*
   *    for the HPDD we need to deal with the truncation
   *     =>  in sampling, if empty topics exist, act as if its 
   *         not truncated, 
   *     =>  once all topics in use ensure normalisation
   */
  if ( ddS.TDt[t]==0 ) 
    /*
     *   will only encounter on new topic;
     *   rest are forced to zero by external logic
     */
    return  (ddP.b0+ddP.a0*ddS.TDTnz) / ((double)ddS.TDT+ddP.b0);
  else if ( ddS.TDTnz==ddN.T ) 
    /*
     *  this cases fudges the situation when *all* topics in use
     *  so we normalise differently;
     *  i.e.    the empty topic at K+1 is inaccessible so we
     *          cannot make it possible
     */
    return ((double)ddS.TDt[t]-ddP.a0)/((double)ddS.TDT-ddN.T*ddP.a0);
  /*
   *   normal situation when empty topics exist
   */
  return  ((double)ddS.TDt[t]-ddP.a0)/((double)ddS.TDT+ddP.b0);	 
}

/*
 *  Base distribution for beta
 *
 *  care must be taken when using this:
 *  in H_PDD case,
 *  if we have zero or multiple cases of TwT[j]==0 
 *  then it wouldn't normalise ... so we arrange things to make it
 *  NB.   no special case of zeroing subsequent null words
 *        as per alphabasetopicprob()
 */
double betabasewordprob(int j) {
  if ( ddP.PYbeta==H_PDP) 
    return ddP.betapr[j];
  if ( ddP.PYbeta==H_HDP ) {
    assert(ddP.aw0==0);
    return ((double)ddS.TwT[j]+ddP.bw0*ddP.betapr[j])/((double)ddS.TWT+ddP.bw0);
  } 
  if ( ddS.TwT[j]==0 ) {
    /*
     *    the GEM remainder is distributed evenly over
     *    all the NULL topics
     */
    assert(ddN.W-ddS.TWTnz>0);
    return  (ddP.bw0+ddP.aw0*(double)ddS.TWTnz)
      /((double)(ddN.W-ddS.TWTnz))/((double)ddS.TWT+ddP.bw0);
  } else if ( ddN.W==ddS.TWTnz ) {
    /*
     *  there are no NULL topics so normalise differently
     */
    return ((double)ddS.TwT[j]-ddP.aw0)/
      ((double)ddS.TWT-ddP.aw0*(double)ddN.W);
  }
  /*
   *   regular GEM estimate
   */
  return((double)ddS.TwT[j]-ddP.aw0)/((double)ddS.TWT+ddP.bw0);
}

/*
 *    prob. the doc table indicator is increased, but not forced
 */
void doctableindicatorprob(int d, int t, int Ttot,
			   double *uone, double *uzero) {
  int nn = ddS.Ndt[d][t];
  int tt = ddS.Tdt[d][t];
  double e1, e0;
  e1 = S_UV(ddC.SX,nn,tt+1);
  if ( tt==1 )
    e0 = nn - ddP.apar;
  else
    e0 = S_U(ddC.SX,nn,tt);
  *uone = e1 * (ddP.bpar+ddP.apar*Ttot) * alphabasetopicprob(t) 
    * (tt+1)/(nn+1);
  *uzero = e0 * (nn-tt+1)/(nn+1);
}

#ifdef NG_SCALESTATS
/*
 *    prob. the NG table indicator is increased, but not forced
 */
void NGtableindicatorprob(int d, int t, int Ttot,
			   double *uone, double *uzero) {
  int nn = ddS.Ndt[d][t];
  int tt = ddS.Tdt[d][t];
  double e1, e0;
  e1 = S_UV(ddC.SX,nn,tt+1);
  if ( tt==1 )
    e0 = nn;
  else
    e0 = S_U(ddC.SX,nn,tt);
  *uone = e1 * 
    ((ddP.ngash + ddS.TDt[t]) /(1/ddP.ngasc + ddS.NGscalestats[t]))
    * (tt+1)/(nn+1) / (ddS.UN[d]+ddP.NGbeta[t]);
  *uzero = e0 * (nn-tt+1)/(nn+1) / (ddS.UN[d]+ddP.NGbeta[t]);
}
#endif

/*
 *    the word table indicator is increased, but not forced
 */
void wordtableindicatorprob(int j, int t, double *uone, double *uzero) {
  int nn = ddS.Nwt[j][t];
  int tt = ddS.Twt[j][t];
  double e1, e0;
  /*
   *   fudge to handle multi-threading  case where constraints violated
   */
  if ( nn>0 && tt==0 ) tt=1;
  if ( tt>nn ) tt = nn;
  if ( nn==0 ) {
    /*
     *    should only happen during multi-threading
     */
    tt = 0;
    e1 = 1;
    e0 = 0;
  } else {
    e1 = S_UV(ddC.SY,nn,tt+1);
    if ( tt==1 )
      e0 = nn - ddP.awpar;
    else
      e0 = S_U(ddC.SY,nn,tt);
  }
  *uone = e1 * (ddP_bwpar(t)+ddP.awpar*ddS.TWt[t]) * betabasewordprob(j) 
    * (tt+1)/(nn+1);
  *uzero = e0 * (nn-tt+1)/(nn+1);
}

/*
 *   probability of topic given document
 *
 *       *zerod - set to zero if a new topic is suggested
 *       *tip - set to prob. indicator would be 1
 *       Ttot - total tables
 */
double topicfact(int d, int t, int Ttot, uint16_t *zerod, float *tip) {
  if ( ddP.theta ) 
    return ddP.theta[d][t];
  if ( ddP.PYalpha ) {
    double p;
    if ( ddP.PYalpha==H_HPDD && ddS.TDt[t]==0 ) {
      /*
       *  special case for HPDD with a topic with 0 occupancy
       *  to handle the introduction of a new topic into a document
       */
      if ( *zerod ) {
	/* want to only do first time  */
	p = (ddP.bpar+ddP.apar*Ttot) *
	  (ddP.b0+ddP.a0*ddS.TDTnz)/(ddP.b0+ddS.TDT);
	*zerod = 0;
      } else {
	/*  subsequent times we set it to zero  */
	p = 0;
      }
      *tip = 1.0;
      return p;
    }
    if ( ddS.Tdt[d][t]==0 ) {
#ifndef NDEBUG
      if ( ddS.Ndt[d][t]>0 ) {
	check_Ndt(d);
	assert(ddS.Ndt[d][t]==0);
      }
#endif
      p = ((double)ddP.bpar+ddP.apar*Ttot) * alphabasetopicprob(t);
      *tip = 1.0;
    } else {
      double uone, uzero;
#ifndef NDEBUG
      if ( ddS.Ndt[d][t]==0 ) {
	check_Ndt(d);
	assert(ddS.Ndt[d][t]>0);
      }
#endif      
      doctableindicatorprob(d, t, Ttot, &uone, &uzero);
      p = uone + uzero;
      *tip = uone/(uone + uzero);
    }
    return p;
  } else if ( ddP.PYalpha==H_NG ) {
#ifdef NG_SCALESTATS
      {
	double p, uone, uzero;
	if ( ddS.Tdt[d][t]==0 ) {
	  /*
	   *   topic empty for this doc
	   */
	  *tip = 1.0;
	  return (ddP.ngash + ddS.TDt[t]) /
	    ( (1/ddP.ngasc + ddS.NGscalestats[t]) *
	      ((double)ddS.UN[d]+ddP.NGbeta[t]) );
	}      
	NGtableindicatorprob(d, t, Ttot, &uone, &uzero);
	p = uone + uzero;
	*tip = uone/(uone + uzero);
	return p;
      }
#else
#ifdef NG_SPARSE
      if ( M_docsparse(d,t)==0 )
	return ((ddP.ngs1+ddS.sparseD[t])/(ddN.DTused-ddS.sparseD[t]+ddP.ngs0))
	  *  ddP.NGalpha[t] / ((double)ddS.UN[d]+ddP.NGbeta[t]);
      else
#endif
	return ((double)ddS.Ndt[d][t]+ddP.NGalpha[t])
	  / ((double)ddS.UN[d]+ddP.NGbeta[t]);
      
#endif
  }
  return ((double)ddS.Ndt[d][t]+ddP.alphapr[t]);
}

/*
 *   only used in estimation
 */
double topicprob(int d, int t, int Ttot) {
  if ( ddP.theta ) 
    return ddP.theta[d][t];
  if ( ddP.PYalpha==H_NG ) {
    /*
     *  this isn't normalised
     */
#ifdef NG_SCALESTATS
    double alphaprt;
    alphaprt = (ddP.ngash+ddS.TDt[t])/(1/ddP.ngasc+ddS.NGscalestats[t]);
    return ((double)ddS.Ndt[d][t]+alphaprt)
      / (ddS.UN[d]+ddP.NGbeta[t]);
#else
#ifdef NG_SPARSE
    if ( ddS.Ndt[d][t]==0 )
      return ((ddP.ngs1+ddS.sparseD[t])/(ddN.DTused+ddP.ngs1+ddP.ngs0))
	*  ddP.NGalpha[t] / ((double)ddS.UN[d]+ddP.NGbeta[t]);
    else
#endif
      return ((double)ddS.Ndt[d][t]+ddP.NGalpha[t])
	/ ((double)ddS.UN[d]+ddP.NGbeta[t]);
#endif
  }
  if ( ddP.PYalpha==H_None ) 
    return ((double)ddS.Ndt[d][t]+ddP.alphapr[t])
      / ((double)ddS.NdT[d]+ddP.alphatot);
  if ( ddP.PYalpha==H_HPDD && ddS.TDt[t]==0 ) {
    /*
     *  special case for HPDD with a topic with 0 occupancy
     *  to handle the introduction of a new topic into a document;
     *  spread probability over all possible NULL cases
     */
    return (ddP.bpar+ddP.apar*Ttot)/((double)ddP.bpar+ddS.NdT[d])
      * (ddP.b0+ddP.a0*ddS.TDTnz)/(ddP.b0+ddS.TDT)
      /(ddN.T-ddS.TDTnz);
  }
  /*
   *    standard PYP result using parent base rate
   */
  if ( ddS.Tdt[d][t]==0 ) {
    return ((double)ddP.bpar+ddP.apar*Ttot) * alphabasetopicprob(t)
      / ((double)ddP.bpar+ddS.NdT[d]);
  } 
  return (ddS.Ndt[d][t] - ddP.apar*ddS.Tdt[d][t]
	  + ((double)ddP.bpar+ddP.apar*Ttot) * alphabasetopicprob(t))
    / ((double)ddP.bpar+ddS.NdT[d]);
}

/*
 *   normally the normaliser term is forgotten for topicfact(),
 *   but we need to add it in when doing the doc PYP
 */
double topicnorm(int d) {
  if ( ddP.PYalpha ) 
    return ((double)ddS.NdT[d]+ddP.bpar);
  return ((double)ddS.NdT[d]+ddP.alphatot);
}

/*
 *   probability of word given topic used for estimation,
 *   which means the topic parts are fixed!
 *   unseen words can still be in the topic, just with
 *   less probability
 */
double wordprob(int j, int t) {
  if ( ddP.phi!=NULL )
    return ddP.phi[t][j];
  if ( ddP.PYbeta ) {
    double pnew = ((double)ddP_bwpar(t)+ddP.awpar*ddS.TWt[t]) * betabasewordprob(j);
    double pold = 0;
    if ( ddS.Nwt[j][t]>0 ) 
      pold = (double)ddS.Nwt[j][t]-ddS.Twt[j][t]*ddP.awpar;
    return (pnew+pold)/((double)ddS.NWt[t]+ddP_bwpar(t));
  }
  return ((double)ddS.Nwt[j][t]+ddP.betapr[j]) / ((double)ddS.NWt[t]+ddP.betatot);
}

/*
 *   post. prob of word given topic used in sampling
 */
double wordfact(int j, int t, float *tip) {
  if ( ddP.phi!=NULL ) {
    assert(ddP.phi);
    return ddP.phi[t][j];
  } 
  if ( ddP.PYbeta ) {
    double p;
    if ( ddS.Twt[j][t]==0 ) {
#ifndef NDEBUG
      if ( ddS.Nwt[j][t]>0 ) {
	yap_message("ddS.Nwt[%d][%d]==%d\n", j, t, ddS.Nwt[j][t]);
	assert(ddS.Nwt[j][t]==0);
      }
#endif
      p = ((double)ddP_bwpar(t)+ddP.awpar*ddS.TWt[t]) * betabasewordprob(j);
      *tip = 1.0;
    } else {
      double uone, uzero;
      wordtableindicatorprob(j, t, &uone, &uzero);
      p = uone + uzero;
      *tip = uone/(uone + uzero);
    }
    return p/((double)ddS.NWt[t]+ddP_bwpar(t));
  }
  return ((double)ddS.Nwt[j][t]+ddP.betapr[j]) / ((double)ddS.NWt[t]+ddP.betatot);
}

/*****************************************************
 *           document level PY
 */

/*
 *   probability of topic from doc level PDP  (ddP.bdk!=NULL version)
 *
 *       t    - topic
 *       (i,mi)    -  word index and corresponding multi version
 *       pK   - input contribution to posterior from adding word with topic
 *       *dip - set to prob. indicator would be 1, if NULL leave
 */
double docfact(D_MiSi_t *dD, int t, int i, int mi, double pK, float *dip) {
  int N = dD->Mi[t], S = dD->Si[t];
  int n, s;
  assert(dip);
  *dip = 1;
  if ( ddP.bdk==NULL ) 
    return pK;
  if ( M_multi(i) ) {
    int mii;
    // assert(mi<ddM.dim_multiind || did==ddN.D-1);
    mii = ddM.multiind[mi]-dD->mi_base;
    assert(mii>=0);
    assert(mii<ddM.MI_max);
    n = dD->Mik[mii][t];
    s = dD->Sik[mii][t];
  } else {
    n = s = 0;
  }  
  if ( s==0 ) {
    return pK * (ddP.bdk[t]+ddP.ad*S)/(ddP.bdk[t]+N); 
  } else {
    double one = pK * (ddP.bdk[t]+ddP.ad*S) * (s+1.0)/(n+1.0);
    double zero = (n-s+1.0)/(n+1.0);
    one *= S_UV(ddC.SD,n,s+1);
    if ( s==1 )
      zero *= n - ddP.ad;
    else
      zero *= S_U(ddC.SD,n,s);    
    *dip = one/(one + zero);
    return (one + zero) /(ddP.bdk[t]+N);
  }
  return 0;
}

/*
 *   counterpart to docfact()
 *   only used in estimation, (ddP.bdk!=NULL version)
 */
double docprob(D_MiSi_t *dD, int t, int i, int mi, double pw) {
  int N = dD->Mi[t], S = dD->Si[t];
  int n, s;
  if ( M_multi(i) ) {
    int mii;
    // assert(mi<ddM.dim_multiind || did==ddN.D-1);
    mii = ddM.multiind[mi]-dD->mi_base;
    assert(mii>=0);
    assert(mii<ddM.MI_max);
    n = dD->Mik[mii][t];
    s = dD->Sik[mii][t];
  } else {
    n = s = 0;
  }  
  if ( s==0 ) {
    return pw * (ddP.bdk[t]+ddP.ad*S)/(ddP.bdk[t]+N); 
  } 
  return (pw * (ddP.bdk[t]+ddP.ad*S) + (n-ddP.ad*s))
    /(ddP.bdk[t] + N);
}
