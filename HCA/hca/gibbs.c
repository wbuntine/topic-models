/*
 * Gibbs sampler.
 * Copyright (C) 2011-2014 Wray Buntine
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
#include "data.h"
#include "stats.h"
#include "check.h"
#include "diag.h"
#include "atomic.h"

/*
 *   remove topic, so update all statistics;
 *
 *   if ddP.bdk active, and indicator not set, only do for docXtopic side
 *   this is indicated by setting wid<0
 *
 *   return non-zero if fails due to constraints
 */
int remove_topic(int i, int did, int wid, int t, int mi, int *Td_, D_MiSi_t *dD) {
  char ud = 0;       /*  indicator for docXtopic's */
  char uw = 0;       /*  indicator for docXwords's */
  /*
   *    this doc contributes data to TM (not just doc PYP)
   */
  // check_Ndt(did);
  if ( ddP.PYalpha && 
       ((ddS.Ndt[did][t]==1) ||
	ddS.Tdt[did][t]>ddS.Ndt[did][t]*rng_unit(rngp) ) )
    ud = 1;
  if ( ddP.PYbeta && ddP.phi==NULL && wid>=0 &&
       ((ddS.Nwt[wid][t]==1) ||
	ddS.Twt[wid][t]>ddS.Nwt[wid][t]*rng_unit(rngp) ) )
    uw = 1;
  if ( (ud==1                    //  no changes to table counts
	&& ddS.Ndt[did][t]>1    //  both n and t zero, so safe
	&& ddS.Tdt[did][t]==1 )    //  t>1 so safe to decrease
       ||
       (uw==1                   //  no changes to table counts
	&& ddS.Nwt[wid][t]>1    // etc.
	&& ddS.Twt[wid][t]==1) ) 
    return 1;
  
  if ( PCTL_BURSTY() && misi_blocked(dD, i, mi, t) )
    return 1;
 
  /*
   *  OK to change this topic since will leave
   *  Ndt[did][t] & Tdt[did][t] constraints OK when we remove
   */

  if ( PCTL_BURSTY() )
    misi_decr(dD, i, mi, t, wid);
  
  /*
   *  do doc X topic updates
   */
  ddS.NdT[did]--;
  assert(ddS.Ndt[did][t]>0);
  ddS.Ndt[did][t]--;
  if ( ud ) {
    /*
     *    subtract affect of table indicator for topic PYP
     */
    unfix_tableidtopic(did,t);
    (*Td_)--;
  }

  /*
   *  do topic X word updates
   */

  if ( ddP.phi==NULL && wid>=0 ) {
    atomic_decr(ddS.NWt[t]);
    // assert(ddS.Nwt[wid][t]>0);
    atomic_decr(ddS.Nwt[wid][t]);
  }
  if ( uw ) {
    /*
     *    subtract affect of table indicator for word PYP
     */
    unfix_tableidword(wid,t);
  }
  return 0;
}

/*
 *   got a new topic, so update all statistics;
 *
 *   if ddP.bdk active, and indicator not set, only do for docXtopic side
 *   this is indicated by setting wid<0
 *
 *   assumes .z[i] fully set
 */
void update_topic(int i, int did, int wid, int t, int mi, int *Td_,
		  D_MiSi_t *dD, float ttip, float wtip, float dtip) {
  ddS.Ndt[did][t]++;
  ddS.NdT[did]++;
  if ( PCTL_BURSTY() ) 
    wid = misi_incr(dD, i, mi, t, wid, dtip);	
  /*
   *   figure out reassigning table id
   */
  if ( ddP.PYalpha && 
       (ddS.Ndt[did][t]==1 || rng_unit(rngp) < ttip) ) {
    (*Td_)++;
    fix_tableidtopic(did, t);
  }
  if ( ddP.phi==NULL && wid>=0 ) {
    int val;
    atomic_incr(ddS.NWt[t]);
    val = atomic_incr(ddS.Nwt[wid][t]);
    if ( ddP.PYbeta ) {
      /*
       *   figure out reassigning table id for word PYP
       */
      if ( val==1 || rng_unit(rngp) < wtip ) {
	/*
	 *  we have a new table for the word matrix
	 */
	fix_tableidword(wid,t);
      }
    }
  }
}

//================
// Gibbs sampler
//================

/********************************
 *   code for LDA 
 *****************************/
double gibbs_lda(/*
		  *  fix==GibbsNone for standard ML training/testing
		  *  fix==GibbsHold for word hold-out testing,
                  *       same as GibbsNone but also handles
		  *       train and test words differently
		  */
		 enum GibbsType fix,
		 int Tmax,   //  to put a choke on growth of T 
		 int did,    //  document index
		 int words,  //  do this many
		 float *p,    //  temp store, at least 4*T
		 D_MiSi_t *dD
		 ) {
  int Td_ = 0;
  int i, wid, t, mi;
  double Z, tot;
  double logdoc = 0;
  int StartWord = ddD.NdTcum[did];
  int EndWord = StartWord + words;
  float *wtip = NULL;
  float *ttip = NULL;
  float *dtip = NULL;

  /*
   *   some of the latent variables are not sampled
   *   are kept in the testing version, uses enum GibbsType
   *      fix = global document setting
   *      fix_doc = settings for word in this doc
   *
   *   NB.   if fix==GibbsNone, then fix_doc==fix
   *         if fix==GibbsHold then fix_doc==GibbsHold or GibbsNone
   */
  enum GibbsType fix_doc = fix;

  if ( PCTL_BURSTY() ) {
    mi = ddM.MI[did];
    assert(ddM.multiind[mi]<ddM.dim_Mi);
    assert(mi<ddM.dim_multiind || did==ddN.D-1);
  }
  if ( ddP.PYalpha )
    Td_ = comp_Td(did);
  /*
   *  assign memory for table indicator probabilities out of p[]
   */
  dtip = p+ddN.T;
  wtip = p+2*ddN.T;
  ttip = p+3*ddN.T;

  for (i=StartWord; i<EndWord; i++) {
    uint16_t zerod=1;      
    
    if ( fix==GibbsHold ) {
      if ( pctl_hold(i) )
	fix_doc = GibbsHold;  //   this word is a hold out
      else
	fix_doc = GibbsNone;
    }
    wid=ddD.w[i]; 
    if ( ddS.TDTnz>=Tmax )
      zerod = 0;
    {
      /*******************
       *   first we remove affects of this word on the stats
       *******************/
      t = Z_t(ddS.z[i]); 
      if ( fix_doc!=GibbsHold  ) {
#ifdef TRACE_WT
	if ( wid==TR_W && t==TR_T )
	  yap_message("remove_topic(w=%d,t=%d,d=%d,l=%d, z=%d, N=%d,T=%d)\n",
		      wid,t,did, i, (int)ddS.z[i],
		      (int)ddS.Nwt[wid][t],(int)ddS.Twt[wid][t]);
#endif
	if ( remove_topic(i, did, 
			  (ddP.bdk==NULL||Z_issetr(ddS.z[i]))?wid:-1,
			  t, mi, &Td_, dD) ) {
	  /*
	   *   not allowed, so no stats altered
	   *   so abandon this word, but still keep diagnostics
	   */    
	  if ( ddG.docode && G_isword(wid) ) {
	    double bdterm = 1;
	    int n = sparsemap_word(wid);
	    if ( PCTL_BURSTY() )
	      bdterm = (ddP.bdk[t]+ddP.ad*dD->Si[t])/(ddP.bdk[t]+dD->Mi[t]);
	    ddG.code[n][t] += bdterm;
	  }
	  if ( ddG.doprob ) {
	    double bdterm = 1;
	    if ( 0 && PCTL_BURSTY() )
	      bdterm = (ddP.bdk[t]+ddP.ad*dD->Si[t])/(ddP.bdk[t]+dD->Mi[t]);
	    if ( did<ddN.DT )
	      ddG.prob[did][t] += bdterm;
	    else
	      ddG.tprob[did-ddN.DT][t] += bdterm;
	  }
	  goto endword;
	}
      }
    }
#ifdef TRACE_WT
    if ( wid==TR_W && t==TR_T)
      yap_message("after remove_topic(w=%d,t=%d,d=%d,l=%d,z=%d,N=%d,T=%d)\n",
		  wid,t,did,i, (int)ddS.z[i],
		  (int)ddS.Nwt[wid][t],(int)ddS.Twt[wid][t]);
#endif
    /***********************
     *    get topic probabilities
     ***********************/
    for (t=0, Z=0, tot=0; t<ddN.T; t++) {
	/*
	 *    (fix_doc==GibbsHold) => 
	 *          doing estimation, not sampling so use prob versions
	 *     (fix_doc!=GibbsHold) => 
	 *          doing sampling so use fact versions
	 */
	double tf = (fix_doc==GibbsHold)?topicprob(did,t,Td_):
	  topicfact(did, t, Td_, &zerod, &ttip[t]);
	if ( tf>0 ) {
	  double wf = (fix_doc==GibbsHold)?wordprob(wid, t):
	    wordfact(wid, t, &wtip[t]);
	  tot += tf;
	  if ( ddP.bdk!=NULL ) 
	    wf = (fix_doc==GibbsHold)?docprob(dD, t, i, mi, wf):
              docfact(dD, t, i, mi, wf, &dtip[t]);
	  Z += p[t] = tf*wf;
	} else
	  p[t] = 0;
    }
    /***********************
     *    do any diagnostics/stats collection;
     *    both these collect prob. of topic, so don't
     *    need anything other than standard p[] calcs
     ***********************/    
    if ( ddG.docode && G_isword(wid) ) {
      double bdterm = 1;
      int n = sparsemap_word(wid);
      for (t=0; t<ddN.T; t++) 
	ddG.code[n][t] += p[t]*bdterm/Z;
    }
    if ( ddG.doprob ) {
      double bdterm = 1;
      if ( did<ddN.DT )
	for (t=0; t<ddN.T; t++) 
	  ddG.prob[did][t] += p[t]*bdterm/Z;
      else
	for (t=0; t<ddN.T; t++) 
	  ddG.tprob[did-ddN.DT][t] += p[t]*bdterm/Z;
    }
    if ( fix!=GibbsHold || fix_doc==GibbsHold )
      logdoc += log(Z/tot);
    yap_infinite(logdoc);
    
    /*******************
     *   now sample t using p[] and install affects of this on the stats;
     *   but note this needs indicator to be set!
     *******************/
    if ( fix_doc!=GibbsHold ) {
      /*
       *  sample and update core stats 
       */
      t = samplet(p, Z, ddN.T,rng_unit(rngp) );
      Z_sett(ddS.z[i],t);
#ifdef TRACE_WT
      if ( wid==TR_W && t==TR_T )
	yap_message("update_topic(w=%d,t=%d,d=%d,l=%d,z=%d,N=%d,T=%d)\n",
		    wid,t,did,i,ddS.z[i],
		    (int)ddS.Nwt[wid][t],(int)ddS.Twt[wid][t]);
#endif
      update_topic(i, did, wid,
		   t, mi, &Td_, dD, ttip[t], wtip[t], dtip[t]);
#ifdef TRACE_WT
      if ( wid==TR_W && t==TR_T )
	yap_message("after update_topic(w=%d,t=%d,d=%d,l=%d,z=%d,N=%d,T=%d)\n",
		    wid,t,did,i,ddS.z[i],
		    (int)ddS.Nwt[wid][t],(int)ddS.Twt[wid][t]);
#endif
    }

    endword:
    if ( PCTL_BURSTY() && M_multi(i) ) {
      mi++;
    }
  }
  return logdoc;
}

