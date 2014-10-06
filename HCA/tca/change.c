/*
 * Changing data structures for statistics
 * Copyright (C) 2010-2014 Wray Buntine  and Jinjing Li
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@monash.edu)
 *         Jinjing Li <jinjingli@gmail.com>
 *
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
#include "pctl.h"
#include "atomic.h"

/*
 *  coount non-zero entries in m_vte[t][]
 */
int nonzero_m_vte(int e, int t) {
  int w;
  int nz=0;
  for (w=0; w<ddN.W; w++ )
    if ( ddS.m_vte[w][t][e]>0 )
      nz++;
  return nz;
}

/*
 *  coount non-zero entries in n_dt[][t]
 */
int nonzero_n_dt(int t) {
  int d;
  int nz=0;
  for (d=0; d<ddN.DT; d++ )
    if ( ddS.n_dt[d][t]>0 )
      nz++;
  return nz;
}

/*
 *   find total t's in document
 */
uint16_t comp_Td(int did) {
  uint16_t Td_ = 0;
  int t;
  assert(ddS.c_dt);
  for (t=0; t<ddN.T; t++) {
    Td_ += ddS.c_dt[did][t];
  }
  return Td_;
}

void unfix_tableidtopic(int d, int t, int ind) { 
  int e = ddD.e[d];
#ifndef H_THREADS
  assert(ddS.c_dt[d][t]>0);
#endif
  if ( ddS.c_dt[d][t]==1 ) {
    ddS.c_dt[d][t] = 0;
    //   should be safe since its associated with ddS.c_dt[d][t]
    if ( atomic_decr(ddS.C_eDt[e][t])>=UINT32_MAX-40 ) {
      /*
       *    shouldn't happen as its a total
       */
      yap_quit("Whoops atomic_decr(ddS.C_eDt[%d][%d])==%u>=UINT32_MAX-40\n",
	       e, t, (unsigned)ddS.C_eDt[e][t]);
    }
    ddS.C_dT[d]--;
    if ( atomic_decr(ddS.C_e[e])>=UINT32_MAX-40 )
      yap_quit("Whoops atomic_decr(ddS.C_e[%d])==%u>=UINT32_MAX-40\n",e,
	       (unsigned)ddS.C_e[e]);
    if ( !ddP.mu ) {
      /*
       *  decrementing up 
       */
      int i;
      int end_e;
      end_e = e-ind;
      for (i = e; i>end_e; i--) {
        if ( atomic_decr(ddS.cp_et[i][t])>=UINT32_MAX-40 ) {
          atomic_incr(ddS.cp_et[i][t]);
          break;
        }       
        if ( atomic_decr(ddS.Cp_e[i])>=UINT32_MAX-40 )
          yap_quit("Whoops atomic_decr(ddS.Cp_e[%d])==%u>=UINT32_MAX-40\n", 
                   i, (unsigned)ddS.Cp_e[i]);
      }
#ifdef MU_CACHE
      mu_side_fact_change(i);
#endif
    }
  } else {
    if ( !ddP.mu ) {
      /*
       *  decrementing down and expect non-zero at the end
       */
      int i;
      int end_e;
      end_e = e-ind;
      for (i=end_e+1; i<=e; i++) {
        if ( atomic_decr(ddS.cp_et[i][t])>=UINT32_MAX-40 ) {
          /*
           *    this may decrement below zero if two independently set ind
           *    so we catch it here
           */
          yap_message("Whoops atomic_decr(ddS.cp_et[%d][%d])==%u>=UINT32_MAX-40\n", i, t,
                      (unsigned) ddS.cp_et[i][t]);
          atomic_incr(ddS.cp_et[i][t]);
        } else {
          if ( atomic_decr(ddS.Cp_e[i])>=UINT32_MAX-40 )
            yap_quit("Whoops atomic_decr(ddS.Cp_e[i])>=UINT32_MAX-40\n");
        }
      }
#ifdef MU_CACHE
      mu_side_fact_change(end_e);
#endif
    }
    ddS.c_dt[d][t]--;
    //   always safe since its associated with ddS.c_dt[d][t]
    if ( atomic_decr(ddS.C_eDt[e][t])>=UINT32_MAX-40 ) {
      /*
       *    shouldn't happen as its a total
       */
      yap_quit("Whoops atomic_decr(ddS.C_eDt[e][t])>=UINT32_MAX-40\n");
    }
    ddS.C_dT[d]--;
    if ( atomic_decr(ddS.C_e[e])>=UINT32_MAX-40 )
      yap_quit("Whoops atomic_decr(ddS.C_e[e])>=UINT32_MAX-40\n");
  }
#if 0
  if ( !ddP.mu ) {
    if ( (ddS.C_eDt[e][t]+((e<ddN.E-1)?ddS.cp_et[e+1][t]:0))<ddS.cp_et[e][t] ) {
      int ii;
      yap_message("\ne=%d, ind=%d, C_eDt=%d+%d, cp_et=%d\n"
                  "ddS.C_eDt[e][t]<ddS.cp_et[e][t]\n", 
                  e, ind, i_eDt, ((e<ddN.E-1)?ddS.cp_et[e+1][t]:0),
                  i_et);
    }
    assert( (ddS.C_eDt[e][t]+((e<ddN.E-1)?ddS.cp_et[e+1][t]:0))
            >= ddS.cp_et[e][t] );
  }
#endif
}

void fix_tableidtopic(int d, int t, int ind) {
  int i;
  int end_e;
  int e = ddD.e[d];
  
  assert(e>=0 && e<ddN.E);

  ddS.C_dT[d]++;
  ddS.c_dt[d][t]++;
  
  atomic_incr(ddS.C_e[e]);
  atomic_incr(ddS.C_eDt[e][t]);
  
  if ( ddP.mu )
    return;
  end_e = e-ind;
  for (i = e; i>end_e; i--) {
    atomic_incr(ddS.Cp_e[i]);
    atomic_incr(ddS.cp_et[i][t]); 
  }
#ifdef MU_CACHE
  mu_side_fact_change(i);
#endif
}

void unfix_tableidword(int e, int w, int t, int ind) {
  int i;
  int lasti=-1;
  assert(e-ind+1>=0);
  assert(ind>0);
  for (i=e; i>e-ind; i--) {
    assert(i>=0);
    if ( (atomic_decr(ddS.s_vte[w][t][i]))>=UINT32_MAX-40 ) {
      /*
       *    this may decrement below zero if two independently set ind
       *    so we catch it here, and quit early
       */
      //WRAY   this tapped!!
      yap_message("Whoops atomic_decr(ddS.s_vte[w][t][i])>=UINT32_MAX-40\n");
      atomic_incr(ddS.s_vte[w][t][i]);
      return;
    } 
    if ( atomic_decr(ddS.S_Vte[t][i])>=UINT32_MAX-40 )
      yap_quit("Whoops atomic_decr(ddS.S_Vte[t][i])>=UINT32_MAX\n");
    lasti = i;
  }    
#ifdef PHI_NORM_CACHE
  phi_norm_change(t,lasti-1);
#endif
#ifdef PHI_CACHE
  phi_unit_change(w,t,lasti-1);
#endif
  if ( lasti==0 ) {
    int val;
#ifndef H_THREADS
    assert(ddS.S_0vT[w]>0);
#endif
    if ( (val=atomic_decr(ddS.S_0vT[w]))>=UINT32_MAX-40 ) {
      yap_message("Whoops atomic_decr(ddS.S_0vT[w])>=UINT32_MAX\n");
      atomic_incr(ddS.S_0vT[w]);
      return;
    }
    if ( atomic_decr(ddS.S_0)>=UINT32_MAX-40 )
      yap_quit("Whoops atomic_decr(ddS.S_0)>=UINT32_MAX %u\n", ddS.S_0);
    if ( val==0 ) {
      if ( atomic_decr(ddS.S_0_nz)>=UINT32_MAX-40 )    
	yap_quit("Whoops atomic_decr(ddS.S_0_nz)>=UINT32_MAX %u\n", ddS.S_0_nz);
    }
  }
}

void fix_tableidword(int e, int w, int t, int ind) { 
  int i;
  int lasti = -1;
  for (i=e; i>e-ind; i--) {
    atomic_incr(ddS.S_Vte[t][i]);
    atomic_incr(ddS.s_vte[w][t][i]);
    lasti = i;
  } 
#ifdef PHI_CACHE
  phi_unit_change(w,t,lasti-1);
#endif
#ifdef PHI_NORM_CACHE
  phi_norm_change(t,lasti-1);
#endif
  if ( lasti==0 ) {
    int val;
    atomic_incr(ddS.S_0);
    val = atomic_incr(ddS.S_0vT[w]);
    if ( val==1 ) {
      atomic_incr(ddS.S_0_nz);    
    }
  }
}

/*
 *    add count to s_vte[][][] stats to make consistent
 *    rippling back if needed
 */
static void add_tableidword(int e, int w, int t) { 
  int laste = -1;
#ifndef H_THREADS
  assert(ddS.s_vte[w][t][e]==0);
#endif
  for (;  e>=0 && ddS.s_vte[w][t][e]==0; e--) {
    int zero = 0;
    //WRAY ???  only increment if zero ... how to do safely
    if ( atomic_incr_val(ddS.s_vte[w][t][e],zero) ) {
      atomic_incr(ddS.S_Vte[t][e]);
      laste = e;
    } else
      break;
  }
  if ( laste==0 ) {
    int val;
    /*   we incremented  ddS.s_vte[w][t][0] */
    atomic_incr(ddS.S_0);
    val = atomic_incr(ddS.S_0vT[w]);
    if ( val==1 ) {
       atomic_incr(ddS.S_0_nz);    
    }
 }
}

/*
 *    remove affects of document from stats
 */
int remove_doc(int d, enum GibbsType fix) {
  int i, t;
  int e = ddD.e[d];
  for (t=0; t<ddN.T; t++) 
    ddS.n_dt[d][t] = 0;
  ddS.N_dT[d] = 0;

  /*
   *    remove topic counts from c & C stats
   */
  for (t=0; t<ddN.T; t++) 
    if ( ddS.c_dt[d][t]>0 ) {
      ddS.C_eDt[e][t] -= ddS.c_dt[d][t];
      /*
       *    remove topic counts from cp stats to make consistent
       *    rippling back if needed
       */
      if ( !ddP.mu ) {
        i = e;
        if ( i==ddN.E-1 ) {
          if ( ddS.C_eDt[i][t]<ddS.cp_et[i][t] ) {
            int diff = ddS.cp_et[i][t] - ddS.C_eDt[i][t];
            ddS.cp_et[i][t] = ddS.C_eDt[i][t];
            ddS.Cp_e[i] -= diff;
          }
          i--;
        }
        for (;  i>=0; i--) {
          int thisc = ddS.C_eDt[i][t] + ddS.cp_et[i+1][t];
          if ( thisc<ddS.cp_et[i][t] ) {
            int diff = ddS.cp_et[i][t]-thisc;
            ddS.cp_et[i][t] = thisc;
            ddS.Cp_e[i] -= diff;
          } else
            break;
        }
      }
      ddS.c_dt[d][t] = 0;
    }      
  ddS.C_e[e] -= ddS.C_dT[d];
  ddS.C_dT[d] = 0;
  
  for (i=ddD.N_dTcum[d]; i<ddD.N_dTcum[d+1]; i++) {
    if ( fix!=GibbsHold || !pctl_hold(i) ) {
      /*
       *   these words are for training
       */
      if ( !PCTL_BURSTY() || Z_issetr(ddS.z[i]) ) {
        int w = ddD.w[i];
        int e1;
        int laste1=-1;
        t = Z_t(ddS.z[i]);
        ddS.M_Vte[t][e]--;
        assert(ddS.m_vte[w][t][e]>0);
        ddS.m_vte[w][t][e]--;
        if ( !ddP.phi ) {
          /*
           *    remove count from s_vte[][][] stats to make consistent
           *    rippling back if needed
           */
          e1 = e;
          if ( e1==ddN.E-1 ) {
            if ( ddS.s_vte[w][t][e1]>ddS.m_vte[w][t][e1] ) {
              ddS.s_vte[w][t][e1]--;
              ddS.S_Vte[t][e1]--;
              laste1 = e1;
            }
            e1--;
          }
          for (;  e1>=0; e1--) {
            int thism = ddS.s_vte[w][t][e1+1]+ddS.m_vte[w][t][e1];
            if ( thism < ddS.s_vte[w][t][e1] ) {
              ddS.s_vte[w][t][e1]--;
              ddS.S_Vte[t][e1]--;
              laste1 = e1;
            } else
              break;
          }
          if ( laste1==0 ) {
            /*   we decremented  ddS.s_vte[w][t][0] */
            ddS.S_0--;
            assert(ddS.S_0vT[w]>0);
            ddS.S_0vT[w]--;
            if ( ddS.S_0vT[w]==0 ) {
              ddS.S_0_nz--;    
            }
          }
        }
      }
    }
  }
  return 0;
}

/*
 *    add affects of document to stats
 *
 *    add minimal indicators/table counts to keep legal
 *
 *    when using multis, reset ddD.Mi[] to be totals
 *    for training words only, ignore test words in hold
 */
int add_doc(int d, enum GibbsType fix) {
  int i, t, w, nd=0;
  int mi = 0;
  int e = ddD.e[d];
  if ( PCTL_BURSTY() ) {
    mi = ddM.MI[d];
  }

  for (i=ddD.N_dTcum[d]; i<ddD.N_dTcum[d+1]; i++) {
    if ( fix!=GibbsHold || pctl_hold(i) ) {
       /*
       *   these words are for perp. calcs
       */     
      nd++;
    }
    if ( fix!=GibbsHold || !pctl_hold(i) ) {
      t = Z_t(ddS.z[i]);
      /*
       *   these words are for training
       */
      ddS.n_dt[d][t]++;
      ddS.N_dT[d]++;  
      if ( !PCTL_BURSTY() || Z_issetr(ddS.z[i]) ) {
        w = ddD.w[i];
        ddS.M_Vte[t][e]++;
        ddS.m_vte[w][t][e]++;
        if ( !ddP.phi && ddS.s_vte[w][t][e]==0 ) 
          add_tableidword(e,w,t);
      }
    }
    if ( PCTL_BURSTY() && M_multi(i) ) mi++;
  }
  /*  initialise ddS.c_dt[d][*]  */ 
  /*
   *   adjust table count stats based on n_dt[d]
   */
  ddS.C_dT[d] = 0;
  for (t=0; t<ddN.T; t++) {
    if ( ddS.n_dt[d][t]>0 ) {
      fix_tableidtopic(d,t,0);
    } else
      ddS.c_dt[d][t] = 0;
  }
  // yap_message("add_doc(%d):  %d\n", d, nd);
  return nd;
}


