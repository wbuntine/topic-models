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
 *  print out the topic topk=10 words. report the PMI score. 
 */
double report_pmi(char *topfile,   /* name of topics file */
		  char *pmifile,  /* name of PMI file */
		  int T,          /* total topics */
		  int W,          /* total words */
		  int E,
		  int topk,
		  double *tp);

/*
 *    computes various quality measures,
 *    sending them both to the log and to the ".par" file
 */
void tca_report(char *resstem, char *stem, int ITER, int procs, 
                enum GibbsType fix, int dopmi) {
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
    if ( dopmi ) {
      char *topfile, *pmifile;
      double coh;
      double *tp;
      tp = dvec(ddN.T);
      topfile=yap_makename(resstem,".top");
      pmifile=yap_makename(stem,".pmi");
      get_probs(tp);
      coh = report_pmi(topfile, pmifile, ddN.T, ddN.W, ddN.E, 10, tp);
      free(topfile);
      free(pmifile);
      if ( fp )
        fprintf(fp, "PMI_%d = %lf\n", ITER, coh);
      free(tp);
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
 *
 *    warm=1  then leave N_dt, N_dT, m_evt and M_eVt alone too
 */
void tca_reset_stats(char *resstem, int restart, int warm) {  
  int e, i, t;
  /*
   *  initialisation *not* done for test docs
   */
  if ( warm==0 ) {
    /*
     *   build the basic stats from z
     */
    memset((void*)ddS.N_dT, 0, sizeof(ddS.N_dT[i])*ddN.D);
    memset((void*)ddS.n_dt[0], 0, sizeof(ddS.n_dt[i][t])*ddN.D*ddN.T);
    memset((void*)ddS.m_evt[0][0], 0, sizeof(ddS.m_evt[e][i][t])*ddN.W*ddN.E*ddN.T);
    memset((void*)ddS.M_eVt[0], 0, sizeof(ddS.M_eVt[e][t])*ddN.E*ddN.T);
    for (i=0; i<ddN.NT; i++) { 
      int d = ddD.d[i];
      t = Z_t(ddS.z[i]);
      ddS.n_dt[d][t]++; 
      ddS.N_dT[d]++; 
      if ( !PCTL_BURSTY() || Z_issetr(ddS.z[i]) ) {
	e = ddD.e[d];
	ddS.m_evt[e][ddD.w[i]][t]++; 
	ddS.M_eVt[e][t]++; 
      }
      assert(ddD.d[i]<ddN.DT);
    }
  }
  memset((void*)ddS.C_dT, 0, sizeof(ddS.C_dT[i])*ddN.D);
  memset((void*)ddS.S_0vT, 0, sizeof(ddS.S_0vT[i])*ddN.W);
  memset((void*)ddS.C_e, 0, sizeof(ddS.C_e[i])*ddN.E);
  memset((void*)ddS.Cp_e, 0, sizeof(ddS.Cp_e[i])*ddN.E);
  memset((void*)ddS.c_dt[0], 0, sizeof(ddS.c_dt[i][t])*ddN.D*ddN.T);
  memset((void*)ddS.S_eVt[0], 0, sizeof(ddS.S_eVt[i][t])*ddN.E*ddN.T);
  memset((void*)ddS.C_eDt[0], 0, sizeof(ddS.C_eDt[i][t])*ddN.E*ddN.T);
  memset((void*)ddS.cp_et[0], 0, sizeof(ddS.cp_et[i][t])*ddN.E*ddN.T);
  memset((void*)ddS.s_evt[0][0], 0, sizeof(ddS.s_evt[e][i][t])*ddN.W*ddN.E*ddN.T);
  ddS.S_0 = 0;
  ddS.S_0_nz = 0;

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
		    (thism<ddS.s_evt[e][i][t]) ) {
	    if ( restart )
	      yap_quit("Inconsistency:  ddS.m_evt[%d][%d][%d]=%d, "
		       "ddS.s_evt[%d][%d][%d]=%d\n",
		       e, i, t, (int)thism,
		     e, i, t, (int)ddS.s_evt[e][i][t]);
	    else {
	      if ( ddS.s_evt[e][i][t]>thism )
		ddS.s_evt[e][i][t] = thism;
	    }
	  }
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
	int n = ddS.n_dt[i][t];
	if ( (n>0) && (ddS.c_dt[i][t]==0) ) 
	  ddS.c_dt[i][t] = 1;
	else if ( ((n>0) ^ (ddS.c_dt[i][t]>0)) ||
		  (ddS.c_dt[i][t]>n)) { 
	  if ( restart ) 
	    yap_quit("Inconsistency:  ddS.n_dt[%d][%d]=%d, ddS.c_dt[%d][%d]=%d\n",
		     i, t, (int)n, i, t, (int)ddS.c_dt[i][t]);
	  else {
	    if ( ddS.c_dt[i][t]>n )
	      ddS.c_dt[i][t] = n;
	  }
	}
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
		  (ddS.cp_et[e][t]>thisc) ) {
	  if ( restart )
	    yap_quit("Inconsistency:  ddS.C_eDt[%d][%d]=%d, "
		     "ddS.cp_et[%d][%d]=%d\n",
		     e, t, (int)thisc, e, t, (int)ddS.cp_et[e][t]);
	  else {
	    if ( ddS.cp_et[e][t]>thisc )
	      ddS.cp_et[e][t] = thisc;
	  }
	}
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


