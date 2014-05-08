/*
 * Various data structures for statistics
 * Copyright (C) 2010-2014 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@monash.edu)
 *
 *     Various data structure read/write/report routines.
 *     Defined in "tca.h"  "stats.h", ...
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
 *    computes various quality measures,
 *    sending them both to the log and to the ".par" file
 */
void tca_report(char *resstem, char *stem, int ITER, int procs, 
                enum GibbsType fix) {
    char *fname;
    double logprob;
    fname = yap_makename(resstem,".par");
    FILE *fp = fopen(fname, "a");
    if ( !fp )
      yap_sysquit("Cannot open output '%s' file:", fname);
    logprob = likelihood();
    fprintf(fp, "logperptrain = %lf\n", -M_LOG2E * logprob/ddN.NT);
    yap_message("log_2(train perp) = %lf\n", -M_LOG2E * logprob/ddN.NT);
    if ( ddN.TEST>0 ) { 
      if ( ddP.mltiter>0 ) {
	char *teststr = fix==GibbsHold?"Hold":"ML";
	logprob = lp_test_ML(fix);
	yap_message("log_2(test perp%s) = %lf\n", teststr, -M_LOG2E * logprob);
	fprintf(fp, "logperp%stest_%d = %lf\n", 
		teststr, ITER, -M_LOG2E * logprob);
      }
    }
    fclose(fp);
    free(fname);
}

/*
 *  allocation of various data statistics, matrices, vectors
 */
void tca_alloc() {
  ddS.z = u16vec(ddN.N);
  if ( !ddS.z )
    yap_quit("Cannot allocate memory for data\n");
  ddS.s_evt = u32tri(ddN.E,ddN.W,ddN.T);
  ddS.S_0vT = u32vec(ddN.W);
  ddS.S_eVt = u32mat(ddN.E,ddN.T);
  ddS.c_dt = u16mat(ddN.D,ddN.T);
  ddS.C_dT = u16vec(ddN.D);
  ddS.C_eDt = u32mat(ddN.E,ddN.T);
  ddS.cp_et = u32mat(ddN.E,ddN.T);
  ddS.C_e = u32vec(ddN.E);
  ddS.Cp_e = u32vec(ddN.E);
  ddS.N_dT = u16vec(ddN.D);
  ddS.n_dt = u16mat(ddN.D,ddN.T);
  ddS.m_evt = u32tri(ddN.E,ddN.W,ddN.T);
  ddS.M_eVt  = u32mat(ddN.E,ddN.T);
}

void tca_free() {
  /*
   *  free
   */
  free(ddS.m_evt[0][0]);  free(ddS.m_evt[0]);  free(ddS.m_evt);
  free(ddS.M_eVt[0]); free(ddS.M_eVt);
  free(ddS.z);
  free(ddS.N_dT);
  free(ddS.n_dt[0]); free(ddS.n_dt);
  free(ddS.S_0vT); 
  free(ddS.S_eVt[0]); free(ddS.S_eVt);
  free(ddS.s_evt[0][0]);  free(ddS.s_evt[0]);  free(ddS.s_evt);
  free(ddS.C_eDt[0]);  free(ddS.C_eDt);
  free(ddS.c_dt[0]);  free(ddS.c_dt);
  free(ddS.C_dT);
  free(ddS.C_e);
  free(ddS.cp_et[0]);  free(ddS.cp_et);
  free(ddS.Cp_e);
}

/*
 *  coount non-zero entries in m_evt[][t]
 */
int nonzero_m_evt(int e, int t) {
  int w;
  int nz=0;
  for (w=0; w<ddN.W; w++ )
    if ( ddS.m_evt[e][w][t]>0 )
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
  
  ddS.c_dt[d][t]--;
  ddS.C_dT[d]--;
#ifndef H_THREADS
  assert(ddS.c_dt[d][t]>=0);
  assert(ddS.c_dt[d][t]>0 || ddS.n_dt[d][t]==0);
#endif
  //   always safe since its associated with ddS.c_dt[d][t]
  atomic_decr(ddS.C_eDt[e][t]);
  atomic_decr(ddS.C_e[e]);  
  {
    int i;
    int end_e;
    end_e = e-ind;
    for (i = e; i>end_e; i--) {
      if ( atomic_decr(ddS.cp_et[i][t])>=UINT32_MAX-40 ) {
	/*
	 *    this may decrement below zero if two independently set ind
	 *    so we catch it here
	 */
	yap_message("Whoops atomic_decr(ddS.cp_et[i][t])>=UINT32_MAX-40\n");
	atomic_incr(ddS.cp_et[i][t]);
	break;
      }           
      atomic_decr(ddS.Cp_e[i]);
    }
  }
}

void fix_tableidtopic(int d, int t, int ind) {
  int i;
  int end_e;
  int e = ddD.e[d];
  
  assert(e>=0 && e<ddN.E);

  ddS.c_dt[d][t]++;
  ddS.C_dT[d]++;
  
  atomic_incr(ddS.C_e[e]);
  atomic_incr(ddS.C_eDt[e][t]);
  
  end_e = e-ind;
  for (i = e; i>end_e; i--) {
    atomic_incr(ddS.cp_et[i][t]); 
    atomic_incr(ddS.Cp_e[i]);
  }
}

void unfix_tableidword(int e, int w, int t, int ind) {
  int i;
  int lasti=-1;
  assert(e-ind+1>=0);
  for (i=e; i>e-ind; i--) {
    assert(i>=0);
    if ( atomic_decr(ddS.s_evt[i][w][t])>=UINT32_MAX-40 ) {
      /*
       *    this may decrement below zero if two independently set ind
       *    so we catch it here
       */
      yap_message("Whoops atomic_decr(ddS.s_evt[i][w][t])>=UINT32_MAX-40\n");
      atomic_incr(ddS.s_evt[i][w][t]);
      break;
    } 
    atomic_decr(ddS.S_eVt[i][t]);
    lasti = i;
  }    
  if ( lasti==0 ) {
    atomic_decr(ddS.S_0);
#ifndef H_THREADS
    assert(ddS.S_0vT[w]>0);
#endif
    atomic_decr(ddS.S_0vT[w]);
    atomic_decr(ddS.S_0_nz);    
  }
}

void fix_tableidword(int e, int w, int t, int ind) { 
  int i;
  int lasti = -1;
  for (i=e; i>e-ind; i--) {
    atomic_incr(ddS.s_evt[i][w][t]);
    atomic_incr(ddS.S_eVt[i][t]);
    lasti = i;
  } 
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
 *    add count to s_evt[][][] stats to make consistent
 *    rippling back if needed
 */
static void add_tableidword(int e, int w, int t) { 
  int laste = -1;
#ifndef H_THREADS
  assert(ddS.s_evt[e][w][t]==0);
#endif
  for (;  e>=0 && ddS.s_evt[e][w][t]==0; e--) {
    //WRAY ???  only increment if zero ... how to do safely
    if ( atomic_incr_zero(ddS.s_evt[e][w][t]) ) {
      atomic_incr(ddS.S_eVt[e][t]);
      laste = e;
    } else
      break;
  }
  if ( laste==0 ) {
    int val;
    /*   we incremented  ddS.s_evt[0][w][t] */
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
        ddS.M_eVt[e][t]--;
        assert(ddS.m_evt[e][w][t]>0);
        ddS.m_evt[e][w][t]--;
	/*
	 *    remove count from s_evt[][][] stats to make consistent
	 *    rippling back if needed
	 */
	e1 = e;
	if ( e1==ddN.E-1 ) {
          if ( ddS.s_evt[e1][w][t]>ddS.m_evt[e1][w][t] ) {
            ddS.s_evt[e1][w][t]--;
            ddS.S_eVt[e1][t]--;
            laste1 = e1;
          }
	  e1--;
	}
	for (;  e1>=0; e1--) {
	  int thism = ddS.s_evt[e1+1][w][t]+ddS.m_evt[e1][w][t];
	  if ( thism < ddS.s_evt[e1][w][t] ) {
	    ddS.s_evt[e1][w][t]--;
	    ddS.S_eVt[e1][t]--;
            laste1 = e1;
	  } else
	    break;
	}
	if ( laste1==0 ) {
	  /*   we decremented  ddS.s_evt[0][w][t] */
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
        ddS.M_eVt[e][t]++;
        ddS.m_evt[e][w][t]++;
        if ( ddS.s_evt[e][w][t]==0 ) 
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

/*
 *    randomise topics for docs in range
 */
void tca_rand_z(int Tinit, int firstdoc, int lastdoc) {
  int i, t;
  int firstword = ddD.N_dTcum[firstdoc];
  int lastword = ddD.N_dTcum[lastdoc];
  for (i=firstword; i<lastword; i++) {
    t = (int)(Tinit*rng_unit(rngp));
    if ( t>=Tinit ) t = Tinit-1;
    //   zero's r as a side effect
    ddS.z[i] = t;
    //    initally, all doc indicators are on
    Z_setr(ddS.z[i]);
  }
#if 0
  if ( PCTL_BURSTY() ) {
    dmi_rand(&ddM, firstdoc, lastdoc);
  }
#endif
}

/*
 *   read topic assignments (z's) from file (both t and r)
 *   for documents docstart to (docend-1)
 */
void tca_read_z(char *resstem, int docstart, int docend) {
  FILE *fr;
  int i, t, r;
  char buf[50];
  char *restartfile = yap_makename(resstem, ".zt");
  fr = fopen(restartfile,"r");
  if ( !fr )
    yap_sysquit("restart file '%s' not read\n", restartfile);
  i = 0;
  docstart = ddD.N_dTcum[docstart];
  docend = ddD.N_dTcum[docend];
  if ( docstart>0 ) {
    for ( ; i<docstart; i++) {
      if ( !fgets(&buf[0], 50, fr) )
	yap_quit("Cannot read %d-th entry from '%s'\n", i, restartfile);
    }
  }
  for (; i<docend; i++) {
    if ( !fgets(&buf[0], 50, fr) )
      yap_quit("Cannot read %d-th entry from '%s'\n", i, restartfile);
    if ( ( PCTL_BURSTY() && sscanf(buf," %d,%d", &t, &r)!=2 )
	 || ( !PCTL_BURSTY() && sscanf(buf," %d", &t)!=1 ) )
      yap_quit("Cannot read %d-th entry from '%s'\n", i, restartfile);
    if ( t>=ddN.T || t<0 )
      yap_quit("Illegal t=%d\n", t);
    ddS.z[i] = t;
    if ( PCTL_BURSTY() && r )
      Z_setr(ddS.z[i]);
  }
  fclose(fr);
  free(restartfile);
}

void tca_write_z(char *resstem) 
{
  char *fname = yap_makename(resstem, ".zt");
  FILE *fp = fopen(fname,"w");
  int i;
  if ( !fp )
    yap_sysquit("Cannot open file '%s' for write\n", fname);
  for (i = 0; i < ddN.N; i++) {
    if ( PCTL_BURSTY() )
      fprintf(fp, "%u,%u\n", (unsigned)Z_t(ddS.z[i]), (unsigned)Z_r(ddS.z[i]));
    else
      fprintf(fp, "%u\n", (unsigned)Z_t(ddS.z[i]));
  }
  if ( ferror(fp) )
    yap_sysquit("Error on writing file '%s' ", fname);
  fclose(fp);
  free(fname);
}

/*
 *    zero everything and rebuild entirely from z[] (both t and r)
 *    but only for training docs
 *    
 *    restart!=0, resstem!+NULL:  set to stem for results to read
 *    OR
 *    restart!=0, resstem==NULL:   recompute from memory
 *    OR
 *    restart==0:  recompute and set table stats to 1 or minimum
 */
void tca_reset_stats(char *resstem, int restart) {  
  int e, i, t;
  /*
   *  initialisation *not* done for test docs
   */
  for (i=0; i<ddN.D; i++) {
    ddS.C_dT[i]=0;
    ddS.N_dT[i]=0;
  }
  for (i=0; i<ddN.W; i++) {
    ddS.S_0vT[i]=0;
  }
  for (e=0; e<ddN.E; e++) {
    ddS.C_e[e] = 0;
    ddS.Cp_e[e] = 0;
  }
  for (t=0; t<ddN.T; t++) {
    for (i=0; i<ddN.D; i++) {
      ddS.n_dt[i][t]=0;
      ddS.c_dt[i][t]=0;
    }
    for (e=0; e<ddN.E; e++) {
      ddS.M_eVt[e][t]=0;
      ddS.S_eVt[e][t]=0;
      for (i=0; i<ddN.W; i++) {
	ddS.m_evt[e][i][t]=0;
	ddS.s_evt[e][i][t]=0;
      }
      ddS.C_eDt[e][t] = 0;
      ddS.cp_et[e][t] = 0;
    }
  }    
  ddS.S_0 = 0;
  ddS.S_0_nz = 0;
  
  for (i=0; i<ddN.NT; i++) { 
    t = Z_t(ddS.z[i]);
    if ( !PCTL_BURSTY() || Z_issetr(ddS.z[i]) ) {
      int d = ddD.d[i];
      e = ddD.e[d];
      ddS.m_evt[e][ddD.w[i]][t]++; 
      ddS.M_eVt[e][t]++; 
      ddS.n_dt[d][t]++; 
      ddS.N_dT[d]++; 
    }
    assert(ddD.d[i]<ddN.DT);
  }

  if ( restart ) {
    if ( resstem ) {
      char *fname = yap_makename(resstem,".sevt");
      read_u32sparsetri(ddN.E,ddN.W,ddN.T,ddS.s_evt,fname);
      free(fname);
    }
    /*  
     *  check and fix consistency,
     *  note have to work backwards because need ddS.s_evt[e+1][i][t]
     *  to be corrected before working on ddS.s_evt[e][i][t]
     */
    for (e=ddN.E-1; e>=0; e--)
      for (i=0; i<ddN.W; i++) {
	for (t=0; t<ddN.T; t++) {
	  uint32_t thism = ddS.m_evt[e][i][t];
	  if ( e<ddN.E-1) 
	    thism += ddS.s_evt[e+1][i][t];
	  if ( (thism>0) && (ddS.s_evt[e][i][t]==0) )
	    ddS.s_evt[e][i][t] = 1;
	  else if ( ((thism>0) ^ (ddS.s_evt[e][i][t]>0)) ||
		    thism<ddS.s_evt[e][i][t] ) 
	    yap_quit("Inconsistency:  ddS.m_evt[%d][%d][%d]=%d, "
		     "ddS.s_evt[%d][%d][%d]=%d\n",
		     e, i, t, (int)thism,
		     e, i, t, (int)ddS.s_evt[e][i][t]);
	}
      }
  } else {
    /*
     *    guess ddS.s_evt at 1;
     *    backwards due to ripple back effect
     */
    for (e=ddN.E-1; e>=0; e--)
      for (i=0; i<ddN.W; i++) {
	for (t=0; t<ddN.T; t++) {
	  uint32_t thism = ddS.m_evt[e][i][t];
	  if ( e<ddN.E-1) 
	    thism += ddS.s_evt[e+1][i][t];
	  if ( thism>0 ) 
	    ddS.s_evt[e][i][t] = 1;
	}
      }
  } 
  /*
   *    compute S_eVt, S_0vT, S_0_nz, S_0
   */
  for (e=0; e<ddN.E; e++)
    for (i=0; i<ddN.W; i++) {
      for (t=0; t<ddN.T; t++) {
	ddS.S_eVt[e][t] += ddS.s_evt[e][i][t];
      }
    }
  for (i=0; i<ddN.W; i++) {
    for (t=0; t<ddN.T; t++) 
      ddS.S_0vT[i] += ddS.s_evt[0][i][t];
    assert(ddS.S_0vT[i]>=0);
    if ( ddS.S_0vT[i]>0 ) {
      ddS.S_0 += ddS.S_0vT[i];
      ddS.S_0_nz++;
    }
  }
  check_m_evt(0);
  /*
   *     alpha part
   */
  if ( restart ) {
    if ( resstem ) {
      char *fname = yap_makename(resstem,".cdt");
      read_u16sparse(ddN.DT,ddN.T,ddS.c_dt,fname);
      free(fname);
    }
    /*  
     *  check consistency
     *    NB.  writing doesn't record c_dt[i][t]==1, makes it 0,
     *          so we have to correct for this!
     */
    for (i=0; i<ddN.DT; i++) {
      for (t=0; t<ddN.T; t++) {
	if ( (ddS.n_dt[i][t]>0) && (ddS.c_dt[i][t]==0) ) 
	  ddS.c_dt[i][t] = 1;
	else if ( ((ddS.n_dt[i][t]>0) ^ (ddS.c_dt[i][t]>0)) ||
		  (ddS.c_dt[i][t]>ddS.n_dt[i][t])) 
	  yap_quit("Inconsistency:  ddS.n_dt[%d][%d]=%d, ddS.c_dt[%d][%d]=%d\n",
		   i, t, (int)ddS.n_dt[i][t], i, t, (int)ddS.c_dt[i][t]);
      }
    }
  } else {
    for (i=0; i<ddN.DT; i++) {
      for (t=0; t<ddN.T; t++) {
	if ( ddS.n_dt[i][t]>0 ) {
	  ddS.c_dt[i][t] = 1;
	}
      }
    }
  }
  for (i=0; i<ddN.DT; i++) {
    e = ddD.e[i];
    for (t=0; t<ddN.T; t++) {
      int tt = ddS.c_dt[i][t];
      if ( tt>0 ) {
	ddS.C_eDt[e][t] += tt;
	ddS.C_dT[i] += tt;
      }
    }
  }
  for (e=0; e<ddN.E; e++)
    for (t=0; t<ddN.T; t++)
      ddS.C_e[e] += ddS.C_eDt[e][t];
  if ( restart ) {
    if ( resstem ) {
      char *fname = yap_makename(resstem,".cpet");
      read_u32sparse(ddN.E,ddN.T,ddS.cp_et,fname);
      free(fname);
    }
    /*  
     *  check and fix consistency
     *    NB.  writing doesn't record cp_et[e][t]==1, makes it 0,
     *          so we have to correct for this!
     *  fix backwards due to rippling back of counts
     */
    for (e=ddN.E-1; e>=0; e--) {
      for (t=0; t<ddN.T; t++) {
	int thisc = ddS.C_eDt[e][t];
	if ( e<ddN.E-1 )
	  thisc += ddS.cp_et[e+1][t];
	if ( (thisc>0) && (ddS.cp_et[e][t]==0) ) 
	  ddS.cp_et[e][t] = 1;
	else if ( ((thisc>0) ^ (ddS.cp_et[e][t]>0)) ||
		  (ddS.cp_et[e][t]>thisc)) 
	  yap_quit("Inconsistency:  ddS.C_eDt[%d][%d]=%d, "
		   "ddS.cp_et[%d][%d]=%d\n",
		   e, t, (int)ddS.C_eDt[e][t], e, t, (int)ddS.cp_et[e][t]);
      }
    }
  } else {
    for (e=ddN.E-1; e>=0; e--) {
      for (t=0; t<ddN.T; t++) {
	int thisc = ddS.C_eDt[e][t];
	if ( e<ddN.E-1 )
	  thisc += ddS.cp_et[e+1][t];
	if ( thisc>0 ) {
	  ddS.cp_et[e][t] = 1;
	}
      }
    }
  }
  for (e=0; e<ddN.E; e++) 
    for (t=0; t<ddN.T; t++) 
      ddS.Cp_e[e] += ddS.cp_et[e][t];
}

void check_n_dt(int d) {
#ifndef NDEBUG
  int t;
  for (t=0; t<ddN.T; t++) {
    if ( ddS.n_dt[d][t]==0 && ddS.c_dt[d][t]>0 )
      assert(ddS.n_dt[d][t]>0 || ddS.c_dt[d][t]==0);
    if ( ddS.n_dt[d][t]>0 && ddS.c_dt[d][t]==0 )
      assert(ddS.n_dt[d][t]==0 || ddS.c_dt[d][t]>0);
  }
#endif
}

void check_m_evt(int e) {
#ifndef NDEBUG
  int v, t;
  for (v=0; v<ddN.W; v++) {
    for (t=0; t<ddN.T; t++) {
      int m = ddS.m_evt[e][v][t] +((e<ddN.E-1)?ddS.s_evt[e+1][v][t]:0);
      int s = ddS.s_evt[e][v][t];
      if ( m==0 && s>0 )
	assert(m>0 || ddS.s_evt[e][v][t]==0);
      if ( m>0 && s==0 )
	assert(m==0 || ddS.s_evt[e][v][t]>0);
    }
  }
  for (t=0; t<ddN.T; t++) {
    int total = 0;
    for (v=0; v<ddN.W; v++)
      total += ddS.m_evt[e][v][t];
    assert(total==ddS.M_eVt[e][t]);
    total = 0;
    for (v=0; v<ddN.W; v++)
      total += ddS.s_evt[e][v][t];
    assert(total==ddS.S_eVt[e][t]);
  }
  if ( e==0 ) {
    int tot=0;
    int totnz=0;
    for (v=0; v<ddN.W; v++) {
      int total = 0;
      for (t=0; t<ddN.T; t++) 
        total += ddS.s_evt[0][v][t];
      assert(total==ddS.S_0vT[v]);
      tot += ddS.S_0vT[v];
      if ( ddS.S_0vT[v]>0 )
        totnz++;
    }
    assert(tot==ddS.S_0);
    assert(totnz==ddS.S_0_nz);        
  }
#endif
}


