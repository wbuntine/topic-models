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
#include "tca.h"
#include "data.h"
#include "stats.h"
#ifdef H_THREADS
#include "atomic.h"
#endif

void checkm_evt(int w, int var);
/*
 *    sample single indicator given counts==n and multiplicity==t
 *    return 0 if OK and set *ind
 *    return 1 if bad constraints
 */
static int sampleindicator(int n, int t, char *ind) {
  *ind = 0;
  if ( ((n==1) ||
	t>n*rng_unit(rngp) ) ) {
    *ind = 1;
    if ( n>1 && t==1 ) 
      return 1;
  }
  return 0;
}

/*
 *    sample doc side indicator for doc and topic
 *    return 0 if OK and set *rd
 *    return 1 if bad constraints
 */
static int resample_doc_side_ind(int d, int t, int *rd) {
  int  e = ddD.e[d];
  char  ud;
  if ( sampleindicator(ddS.n_dt[d][t],ddS.c_dt[d][t],&ud) )
    // cant sample, abondon
    return 1;
  if ( ud ) {
    int ee;
    //sample rd from 0 ...e+1 (e+1 means going back to prior)
    *rd = 0;
    for (ee=e; ee>=0; ee--) {
      int sec_part=0;
      if (ee<ddN.E-1) 
	sec_part=ddS.cp_et[ee+1][t];
      if (sampleindicator(ddS.C_eDt[ee][t]+sec_part,ddS.cp_et[ee][t],&ud)) 
	return 1;
      if (ud) {
	(*rd)++;
      } else {
	break;
      }
    }
  } else
    *rd = -1;
  return 0;
}

/*
 *    sample word side indicator for epoch, word and topic
 *    return 0 if OK and set *rw
 *    return 1 if bad constraints
 */
static int resample_word_side_ind(int e,  int v, int t, int *rw) {
  int ee;
  char uw;
  // sample rw from 0 ...e+1 (e+1 means going back to prior)
  *rw = 0;
  for (ee=e;ee>=0;ee--) {
    int sec_part=0;
    if ( ee<ddN.E-1 ) 
      sec_part=ddS.s_evt[ee+1][v][t];
    if (sampleindicator(ddS.m_evt[ee][v][t]+sec_part,
			ddS.s_evt[ee][v][t],&uw)) 
      return 1;
    if (uw) {
      (*rw)++;
    } else {
      break;
    }
  }
  return 0;
}
/*
 *   remove topic, so update all statistics;
 *   return non-zero if fails due to constraints
 *
 *   JINJING:   you need to fix this
 */
int remove_topic(int i, int did, int wid, int t, int mi) {
  int  rw;           /*  indicator for docXwords's */
  int  rd;           /*  indicator for docXtopic's */
  int  e = ddD.e[did];
  /*
   *    this doc contributes data to TM (not just doc PYP)
   */
  // check_n_dt(did);

  if ( resample_doc_side_ind(did, t, &rd) )
    return 1;
  assert(rd>=-1 && rd<=e+1);

  if (wid>=0 && resample_word_side_ind(e, wid, t, &rw) )
    return 1;
  assert(wid<0 || (rw>=0 && rw<=e+1));
  
  /*
   *  OK to change this topic 
   */
  if ( PCTL_BURSTY() && M_multi(i)) {
    assert(ddD.multiind[mi]<ddD.dim_Mi);
    assert(mi<ddD.dim_multiind);
    ddS.Si[ddD.multiind[mi]]--; 
  }  
  /*
   *    remove from doc stats
   */
  ddS.N_dT[did]--;
  assert(ddS.n_dt[did][t]>0);
  ddS.n_dt[did][t]--;
  if ( rd>=0 ) {
    /*
     *    subtract affect of table indicator for topic PYP
     */
    unfix_tableidtopic(did,t,rd);
  }
  /*
   *  remove from topicXword stats
   */
  if ( wid>=0 ) {
    ddS.M_eVt[e][t]--;
    if ( ddS.m_evt[e][wid][t]==0 )
      yap_message("Decrementing ddS.m_evt[e][%d][%d]\n", e, wid, t);
    assert(ddS.m_evt[e][wid][t]>0);
    ddS.m_evt[e][wid][t]--;
    if ( rw>=1 ) {
      /*
       *    subtract affect of table indicator for word PYP
       */
      unfix_tableidword(e,wid,t,rw);
    }
  }
  return 0;
}

/*
 *   got a new topic, so update all statistics;
 *
 */
void update_topic(int i, int did, int wid, int t, int mi) {
  int e = ddD.e[did];
  int rd;
  int rw;
  if ( PCTL_BURSTY() && M_multi(i) ) {
    ddS.Si[ddD.multiind[mi]]++;
    if ( ddS.Si[ddD.multiind[mi]]>ddD.Mi[ddD.multiind[mi]])
      assert(ddS.Si[ddD.multiind[mi]]<=ddD.Mi[ddD.multiind[mi]]);
  }
  /*
   *   fix up doc side
   */
  rd = doc_side_ind(did, t);
  ddS.n_dt[did][t]++;
  ddS.N_dT[did]++;
  if ( rd>=0 ) {
    fix_tableidtopic(did, t, rd);
  }
  /*
   *   fix up topicXword side
   */
  rw = word_side_ind(e, wid, t);
  ddS.M_eVt[e][t]++;
  ddS.m_evt[e][wid][t]++;
  if ( rw>0 ) {
    /*
     *  we have a new table for the word matrix
     */
    fix_tableidword(e,wid,t,rw);
  }
#ifdef IND_STATS
  if ( did < ddN.DT ) {
    if ( rd>=0 && rd<=e )
      ddP.doc_ind_stats[t][e][rd]++;
    if ( rw<=e )
      ddP.word_ind_stats[t][e][rw]++;
  } 
#endif
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
     int did,    //  document index
     int words,  //  do this many
     float *p,    //  temp store
     D_MiSi_t *dD
     ) {
  int i, wid, t, mi;
  int e;
  double Z, tot;
  double logdoc = 0;
  int StartWord = ddD.N_dTcum[did];
  int EndWord = StartWord + words;

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
  }
  e = ddD.e[did];

  for (i=StartWord; i<EndWord; i++) {
    if ( fix==GibbsHold ) {
      if ( pctl_hold(i) )
	fix_doc = GibbsHold;  //   this word is a hold out
      else
	fix_doc = GibbsNone;
    }
    // check_m_evt(e);
    wid=ddD.w[i]; 
    /*******************
     *   first we remove affects of this word on the stats
     *******************/
    t = Z_t(ddS.z[i]); 
    if ( fix_doc!=GibbsHold && 
	 ( !PCTL_BURSTY() || Z_issetr(ddS.z[i])) ) {
#ifdef TRACE_WT
      if ( wid==TR_W && t==TR_T )
	yap_message("remove_topic(w=%d,t=%d,d=%d,l=%d, z=%d, N=%d,T=%d)\n",
		    wid,t,did, i, (int)ddS.z[i],
		    (int)ddS.m_evt[wid][t],(int)ddS.s_evt[wid][t]);
#endif
      if ( remove_topic(i, did, Z_issetr(ddS.z[i])?wid:-1, t, mi) ) {
	goto endword;
      }
    }
#ifdef TRACE_WT
    if ( wid==TR_W && t==TR_T)
      yap_message("after remove_topic(w=%d,t=%d,d=%d,l=%d,z=%d,N=%d,T=%d)\n",
		  wid,t,did,i, (int)ddS.z[i],
		  (int)ddS.m_evt[e][wid][t],(int)ddS.s_evt[e][wid][t]);
#endif
    /***********************
     *    get topic probabilities
     ***********************/
    // check_m_evt(e);
     for (t=0, Z=0, tot=0; t<ddN.T; t++) {
      if ( fix_doc==GibbsHold ) {
        /*
         *   doing estimation, not sampling so use prob versions
         */
        double tf = doc_side_prob(did,t);
        if ( tf>0 ) {
          tot += tf;
          p[t] = tf * word_side_prob(e, wid, t);
          Z += p[t];
	} else
	  p[t] = 0;
      } else {
	/*
	 *   doing sampling so use fact versions
	 */
	double tf = doc_side_fact(did,t);
	if ( tf>0 ) {
	  tot += tf;
	  /*
	   *    for hold out words need to estimate word probs,
	   *    not sample t's
	   */
	  p[t] = tf * word_side_fact(e, wid, t);
	  Z += p[t];
	} else
	  p[t] = 0;
      }
    }
    if ( PCTL_BURSTY() ) {
      /***********************
       *    handle the front end doc PYP which adjusts word probabilities
       *    and log probs plus sets the indicator with Z_Xsetr()
       ***********************/    
      float dip;
      double new_Z = 0;
      if ( M_multi(i) ) {
	int mii = ddD.multiind[mi];
	assert(mi<ddD.dim_multiind);
	if ( fix_doc!=GibbsHold ) {
	  if ( ddS.Si[mii]>=ddD.Mi[mii] ) 
	    assert(ddS.Si[mii]<ddD.Mi[mii]);
	  new_Z = docfact(did, ddD.Mi[mii]-1, ddS.Si[mii], 
			  Z/topicnorm(did), &dip);
	  if ( dip==1 || dip > rng_unit(rngp) )
	    Z_setr(ddS.z[i]);
	  else
	    Z_unsetr(ddS.z[i]);
	} else {
	  /*
	   *    ddD.Mi[mi], ddS.Si[mi] haven't been decremented because
	   *    this is test data and is not part of training
	   */
	  new_Z = docprob(did, ddD.Mi[mii], ddS.Si[mii], Z/tot);
	}
      } else {
	/*
	 *   these words only exist once in entire document,
	 *   and its in the test set so count is 0 in other
	 */
	if ( fix_doc!=GibbsHold ) {
	  new_Z = docfact(did, 0, 0, Z/topicnorm(did), &dip);
	  Z_setr(ddS.z[i]);
	} else {
	  new_Z = docprob(did, 0, 0, Z/tot);
	}
      }
      if ( fix!=GibbsHold || fix_doc==GibbsHold )
	logdoc += log(new_Z);
    } else
      if ( fix!=GibbsHold || fix_doc==GibbsHold )
	logdoc += log(Z/tot);
    yap_infinite(logdoc);
    
    /*******************
     *   now sample t using p[] and install affects of this on the stats;
     *   but note this needs indicator to be set!
     *******************/
    if ( fix_doc!=GibbsHold && (!PCTL_BURSTY() || Z_issetr(ddS.z[i])) ) {
      /*
       *  sample and update core stats 
       */
      t = samplet(p, Z, ddN.T, rng_unit(rngp));
      Z_sett(ddS.z[i],t);
#ifdef TRACE_WT
      if ( wid==TR_W && t==TR_T )
  yap_message("update_topic(w=%d,t=%d,d=%d,l=%d,z=%d,N=%d,T=%d)\n",
        wid,t,did,i,ddS.z[i],
        (int)ddS.m_evt[e][wid][t],(int)ddS.s_evt[e][wid][t]);
#endif
      update_topic(i, did, wid, t, mi);
#ifdef TRACE_WT
      if ( wid==TR_W && t==TR_T )
  yap_message("after update_topic(w=%d,t=%d,d=%d,l=%d,z=%d,N=%d,T=%d)\n",
        wid,t,did,i,ddS.z[i],
        (int)ddS.m_evt[e][wid][t],(int)ddS.s_evt[e][wid][t]);
#endif
    }

    endword:
    if ( PCTL_BURSTY() && M_multi(i) ) {
      mi++;
    }
  }
  return logdoc;
}

