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
#include "atomic.h"

void checkm_evt(int w, int var);
/*
 *    sample single indicator given counts==n and multiplicity==t
 *    return 0 if OK and set *ind
 *    return 1 if bad constraints
 */
static int sampleindicator(int n, int t, char *ind) {
  if ( n==0 ) {
    /*    
     *   shouldn't happen, but may due to inconsistency;
     *   since its sampling up, its zero but should have been
     *   one,, so assume its one
     */
    *ind = 1;
    return 0;
  }
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
 *
 *    *rd is number of epochs back to take count
 *        0==none,  e+1==all the way to root
 *    max of #rd ddP.back
 */
static int resample_doc_side_ind(int d, int t, int *rd) {
  int  e = ddD.e[d];
  char  ud;
  if ( sampleindicator(ddS.n_dt[d][t],ddS.c_dt[d][t],&ud) )
    // cant sample, abandon
    return 1;
  if ( ud ) {
    int ee;
    //sample rd from 0 ...e+1 (e+1 means going back to prior)
    *rd = 0;
    for (ee=e; ee>=0; ee--) {
      int sec_part = 0;
      if (ee<ddN.E-1) 
	sec_part = ddS.cp_et[ee+1][t];
      if (sampleindicator(ddS.C_eDt[ee][t]+sec_part,ddS.cp_et[ee][t],&ud)) 
	return 1;
      // equality condition means decrement is forced back regardless
      if ( (ee>e-ddP.back || ddS.C_eDt[ee][t]+sec_part<=ddS.cp_et[ee][t])
	   && ud) {
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
 *
 *    *rw==-1  =>  not table up to alpha
 *    else     =>  number of epochs in alpha table goes back
 *                 0==none,  e+1==all the way to root)
 *    max of #rw si ddP.back
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
    // equality condition means decrement is forced back regardless
    if ( (ee>e-ddP.back || ddS.m_evt[ee][v][t]+sec_part<=ddS.s_evt[ee][v][t])
	&& uw) {
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
int remove_topic(int i, int did, int wid, int t, int mi, D_MiSi_t *dD) {
  int  rw = 0;           /*  indicator for docXwords's */
  int  rd;           /*  indicator for docXtopic's */
  int  e = ddD.e[did];
  /*
   *    this doc contributes data to TM (not just doc PYP)
   */
  // check_n_dt(did);

  if ( PCTL_BURSTY() && misi_blocked(dD, i, mi, t) )
    return 1;

  if ( resample_doc_side_ind(did, t, &rd) )
    return 1;
  assert(rd>=-1 && rd<=e+1);

  if (wid>=0 && resample_word_side_ind(e, wid, t, &rw) )
    return 1;
  assert(wid<0 || (rw>=0 && rw<=e+1));
  
  /*
   *  OK to change this topic.
   *  NB.  all misi stats are local document so unaffected
   *       by other threads
   */
  if ( PCTL_BURSTY() )
    misi_decr(dD, i, mi, t, wid);

  /*
   *    remove from doc stats ...
   */
  if ( rd>=0 ) {
    /* 
     *   subtract affect of table indicator for topic PYP;
     *   when multicore, need to maintain consistency of constraints
     */
    int n = ddS.n_dt[did][t];
    assert(n>0);
    if ( n==1 ) {
      /*  chains going to zero, work up  */
      // ??????????  wrong ... its ddS.c_dt[d][t] we should worry about!
      ddS.n_dt[did][t] = 0;
      ddS.N_dT[did]--;
      unfix_tableidtopic(did,t,rd);
    } else {
     /*  chains not going to zero, work down  */
      unfix_tableidtopic(did,t,rd);
      ddS.n_dt[did][t] = n-1;
      ddS.N_dT[did]--;
    }
  } else {
    ddS.n_dt[did][t]--;
    ddS.N_dT[did]--;
  }
  
  /*
   *  remove from topicXword stats;
   *  note wid<0 means burstiness caught case
   */
  if ( wid>=0 ) {
    if ( rw>=1 ) {
      int decr = 0; 
      /*
       *  need to safely catch with atomic ops if it goes to zero
       */
      if ( e==ddN.E-1 || ddS.s_evt[e+1][wid][t]==0 ) {
        int one = 1;
        /*  no data inherited from later epoch */
        if ( atomic_decr_val(ddS.m_evt[e][wid][t],one) ) {
          /*   m_evt[e][wid][t] was decremented to zero */
          atomic_decr(ddS.M_eVt[e][t]);
          decr = 1;
        }
      }
      /*
       *    subtract affect of table indicator for word PYP
       */
      unfix_tableidword(e,wid,t,rw);
      if ( decr==0 ) {
        /*  not decremented to zero above */
#ifndef H_THREADS
        assert(ddS.m_evt[e][wid][t]>0);
#endif
        /*  safe because the count is uniquely associated with this word */
        atomic_decr(ddS.m_evt[e][wid][t]);
        atomic_decr(ddS.M_eVt[e][t]);
      }
    } else {
      /*  safe because the count is uniquely associated with this word */
      atomic_decr(ddS.m_evt[e][wid][t]);
      atomic_decr(ddS.M_eVt[e][t]);
    }
  }
  return 0;
}

/*
 *   got a new topic, so update all statistics;
 *
 */
void update_topic(int i, int did, int wid, int t, int mi, 
                  float dtip, D_MiSi_t *dD) {
  int e = ddD.e[did];
  int rd;
  int rw=0;
  /*
   *   fix up doc side
   */
  rd = doc_side_ind(did, t);
  assert(rd>=-1 && rd<=e+1);
  
  if ( PCTL_BURSTY() ) 
    // set wid negative if the word is bursty
    wid = misi_incr(dD, i, mi, t, wid, dtip);

  if ( rd>=0 ) {
    /*
     *    increment affect of table indicator for topic PYP;
     *    when multicore, need to maintain consistency of constraints
     */
    int n = ddS.n_dt[did][t];
    if ( 0 && n==0 ) {
      /*  ripple to one, so go down chain */
      fix_tableidtopic(did, t, rd);
      ddS.N_dT[did]++;
      ddS.n_dt[did][t] = 1;
    } else {
      /*  increment up chain */
      ddS.N_dT[did]++;
      ddS.n_dt[did][t] = n+1;
      fix_tableidtopic(did, t, rd);
    }
  } else {
    ddS.N_dT[did]++;
    ddS.n_dt[did][t]++;
  }

  /*
   *   fix up topicXword side
   */
  if ( wid>=0 ) {
    rw = word_side_ind(e, wid, t);
    assert(rw>=0 && rw<=e+1);
    if ( rw>0 ) {
      /*
       *    subtract affect of table indicator for word PYP;
       *    when multicore, need to maintain consistency of constraints
       */
      int zero = 0;
      int incr = 0;
      if ( atomic_incr_val(ddS.m_evt[e][wid][t],zero) ) {
        if ( 1 ) {
          atomic_incr(ddS.M_eVt[e][t]);
          incr = 1;
        } else {
          atomic_decr(ddS.m_evt[e][wid][t]);
        }
        /*  cancel increment and do afterwards */
      } else {
        atomic_incr(ddS.M_eVt[e][t]);
        atomic_incr(ddS.m_evt[e][wid][t]);
        incr = 1;
      }
      /*
       *  we have a new table for the word matrix
       */
      fix_tableidword(e,wid,t,rw);
      if ( !incr ) {
        atomic_incr(ddS.M_eVt[e][t]);
        atomic_incr(ddS.m_evt[e][wid][t]);
      }
    } else {
      /*  simple case, order doesn't matter  */
      atomic_incr(ddS.M_eVt[e][t]);
      atomic_incr(ddS.m_evt[e][wid][t]);
    }
  }
#ifdef IND_STATS
  if ( did < ddN.DT ) {
    if ( rd>=0 && rd<=e )
      atomic_incr(ddP.doc_ind_stats[t][e][rd]);
    if ( wid>=0 && rw<=e )
      atomic_incr(ddP.word_ind_stats[t][e][rw]);
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
  int i, wid, t, mi=0;
  int e;
  double Z, tot;
  double logdoc = 0;
  int logdocinf = 0;
  int StartWord = ddD.N_dTcum[did];
  int EndWord = StartWord + words;
  float dtip[ddN.T];

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
    if ( fix_doc!=GibbsHold ) {
#ifdef TRACE_WT
      if ( wid==TR_W && t==TR_T )
	yap_message("remove_topic(w=%d,t=%d,d=%d,l=%d, z=%d, N=%d,T=%d)\n",
		    wid, t,did, i, (int)ddS.z[i],
		    (int)ddS.m_evt[wid][t],(int)ddS.s_evt[wid][t]);
#endif
      if ( remove_topic(i, did, (!PCTL_BURSTY()||Z_issetr(ddS.z[i]))?wid:-1, 
                        t, mi, dD) ) {
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
      /*
       *   (fix_doc==GibbsHold) =>
       *       doing estimation, not sampling so use prob versions
       *    else
       *        doing sampling so use fact versions
       */
      double tf = (fix_doc==GibbsHold)?doc_side_prob(did,t):
        doc_side_fact(did,t);
      if ( tf>0 ) {
        double wf = (fix_doc==GibbsHold)?word_side_prob(e, wid, t):
          word_side_fact(e, wid, t);
        tot += tf;
        if ( PCTL_BURSTY() ) 
          wf = (fix_doc==GibbsHold)?docprob(dD, t, i, mi, wf):
            docfact(dD, t, i, mi, wf, &dtip[t]);
        Z += p[t] = tf * wf;
      } else
        p[t] = 0;
    }
    if ( fix!=GibbsHold || fix_doc==GibbsHold )
      logdoc += log(Z/tot);
    if ( logdocinf==0 ) 
      if ( !finite(logdoc) ) {
	logdocinf++;
	yap_infinite(logdoc);
      }

    /*******************
     *   now sample t using p[] and install affects of this on the stats;
     *   but note this needs indicator to be set!
     *******************/
    if ( fix_doc!=GibbsHold ) {
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
      update_topic(i, did, wid, t, mi, dtip[t], dD);
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

