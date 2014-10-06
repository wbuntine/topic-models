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
#include "pmi.h"

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
      coh = report_pmi(topfile, pmifile, NULL, ddN.T, ddN.W, ddN.E, 10, tp, NULL);
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
  if ( !ddP.phi ) {
    ddS.s_vte = u32tri(ddN.W,ddN.T,ddN.E);
    ddS.S_0vT = u32vec(ddN.W);
    ddS.S_Vte = u32mat(ddN.T,ddN.E);
  } else {
    ddS.s_vte = NULL;
    ddS.S_0vT = NULL;
    ddS.S_Vte = NULL;
  }
  if ( !ddP.mu) {
    ddS.cp_et = u32mat(ddN.E,ddN.T);
    ddS.Cp_e = u32vec(ddN.E);
  } else {
    ddS.cp_et = NULL;
    ddS.Cp_e = NULL;
  }
  ddS.c_dt = u16mat(ddN.D,ddN.T);
  ddS.C_dT = u16vec(ddN.D);
  ddS.C_e = u32vec(ddN.E);
  ddS.C_eDt = u32mat(ddN.E,ddN.T);
  ddS.N_dT = u16vec(ddN.D);
  ddS.n_dt = u16mat(ddN.D,ddN.T);
  ddS.m_vte = u32tri(ddN.W,ddN.T,ddN.E);
  ddS.M_Vte  = u32mat(ddN.T,ddN.E);
#ifdef MU_CACHE
  mu_side_fact_init();
#endif
#ifdef PHI_CACHE
  phi_cache_init();
#endif
}

void tca_free() {
  /*
   *  free
   */
  free(ddS.m_vte[0][0]);  free(ddS.m_vte[0]);  free(ddS.m_vte);
  free(ddS.M_Vte[0]); free(ddS.M_Vte);
  free(ddS.z);
  free(ddS.N_dT);
  free(ddS.n_dt[0]); free(ddS.n_dt);
  if ( ddS.s_vte ) {
    free(ddS.S_0vT); 
    free(ddS.S_Vte[0]); free(ddS.S_Vte);
    free(ddS.s_vte[0][0]);  free(ddS.s_vte[0]);  free(ddS.s_vte);
  }
  free(ddS.C_eDt[0]);  free(ddS.C_eDt);
  free(ddS.c_dt[0]);  free(ddS.c_dt);
  free(ddS.C_dT);
  free(ddS.C_e);
  if ( ddS.cp_et ) {
    free(ddS.cp_et[0]);  free(ddS.cp_et);
    free(ddS.Cp_e);
  }
#ifdef MU_CACHE
  mu_side_fact_free();
#endif
#ifdef PHI_CACHE
  phi_cache_free();
#endif
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
  //  yap_message("READING %ld to %ld\n", docstart, docend);

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
      yap_quit("Cannot scan %d-th entry from '%s'\n", i, restartfile);
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
 *    warm=1  then leave N_dt, N_dT, m_vte and M_Vte alone too
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
    memset((void*)ddS.c_dt[0], 0, sizeof(ddS.c_dt[i][t])*ddN.D*ddN.T);
    if ( !ddP.mu ) 
      memset((void*)ddS.cp_et[0], 0, sizeof(ddS.cp_et[i][t])*ddN.E*ddN.T);
    memset((void*)ddS.m_vte[0][0], 0, sizeof(ddS.m_vte[i][t][e])*ddN.W*ddN.E*ddN.T);
    memset((void*)ddS.M_Vte[0], 0, sizeof(ddS.M_Vte[t][e])*ddN.E*ddN.T);
    if ( !ddP.phi ) 
      memset((void*)ddS.s_vte[0][0], 0, sizeof(ddS.s_vte[i][t][e])*ddN.W*ddN.E*ddN.T);
    for (i=0; i<ddN.NT; i++) { 
      int d = ddD.d[i];
      t = Z_t(ddS.z[i]);
      ddS.n_dt[d][t]++; 
      ddS.N_dT[d]++; 
      if ( !PCTL_BURSTY() || Z_issetr(ddS.z[i]) ) {
	e = ddD.e[d];
	ddS.m_vte[ddD.w[i]][t][e]++; 
	ddS.M_Vte[t][e]++; 
      }
      assert(ddD.d[i]<ddN.DT);
    }
  }
  memset((void*)ddS.C_dT, 0, sizeof(ddS.C_dT[i])*ddN.D);
  memset((void*)ddS.C_eDt[0], 0, sizeof(ddS.C_eDt[i][t])*ddN.E*ddN.T);
  memset((void*)ddS.C_e, 0, sizeof(ddS.C_e[i])*ddN.E);
  if ( !ddP.mu ) 
    memset((void*)ddS.Cp_e, 0, sizeof(ddS.Cp_e[i])*ddN.E);
  if ( !ddP.phi ) {
    memset((void*)ddS.S_Vte[0], 0, sizeof(ddS.S_Vte[t][i])*ddN.E*ddN.T);
    memset((void*)ddS.S_0vT, 0, sizeof(ddS.S_0vT[i])*ddN.W);
  }
  ddS.S_0 = 0;
  ddS.S_0_nz = 0;

  if ( ddP.phi ) {
    ;  // do nothing!
  } else if ( restart ) {
    int inconsistency = 0;
    if ( resstem ) {
      char *fname = yap_makename(resstem,".sevt");
      read_u32sparsetri(ddN.W,ddN.T,ddN.E,ddS.s_vte,fname);
      free(fname);
    }
    /*  
     *  check and fix consistency,
     *  note have to work backwards because need ddS.s_vte[i][t][e+1]
     *  to be corrected before working on ddS.s_vte[i][t][e]
     */
    for (e=ddN.E-1; e>=0; e--)
      for (i=0; i<ddN.W; i++) {
	for (t=0; t<ddN.T; t++) {
	  uint32_t thism = ddS.m_vte[i][t][e];
	  if ( e<ddN.E-1) 
	    thism += ddS.s_vte[i][t][e+1];
	  if ( (thism>0) && (ddS.s_vte[i][t][e]==0) )
	    ddS.s_vte[i][t][e] = 1;
	  else if ( ((thism>0) ^ (ddS.s_vte[i][t][e]>0)) ||
		    (thism<ddS.s_vte[i][t][e]) ) {
	    if ( restart ) {
	      inconsistency ++;
	      if ( inconsistency<10)
		yap_message("Inconsistency:  ddS.m_vte[%d][%d][%d]=%d, "
			    "ddS.s_vte[%d][%d][%d]=%d\n",
			    e, i, t, (int)thism,
			    e, i, t, (int)ddS.s_vte[i][t][e]);
	    } else {
	      if ( ddS.s_vte[i][t][e]>thism )
		ddS.s_vte[i][t][e] = thism;
	    }
	  }
	}
      }
    if (  inconsistency )
      yap_quit("Inconsistencies = %d\n", inconsistency);
  } else {
    /*
     *    guess ddS.s_vte at 1;
     *    backwards due to ripple back effect
     */
    for (e=ddN.E-1; e>=0; e--)
      for (i=0; i<ddN.W; i++) {
	for (t=0; t<ddN.T; t++) {
	  uint32_t thism = ddS.m_vte[i][t][e];
	  if ( e<ddN.E-1) 
	    thism += ddS.s_vte[i][t][e+1];
	  if ( thism>0 ) 
	    ddS.s_vte[i][t][e] = 1;
	}
      }
  } 
  if ( !ddP.phi ) {
    /*
     *    compute S_Vte, S_0vT, S_0_nz, S_0
     */
    for (e=0; e<ddN.E; e++)
      for (i=0; i<ddN.W; i++) {
	for (t=0; t<ddN.T; t++) {
	  ddS.S_Vte[t][e] += ddS.s_vte[i][t][e];
	}
      }
    for (i=0; i<ddN.W; i++) {
      for (t=0; t<ddN.T; t++) 
	ddS.S_0vT[i] += ddS.s_vte[i][t][0];
      assert(ddS.S_0vT[i]>=0);
      if ( ddS.S_0vT[i]>0 ) {
	ddS.S_0 += ddS.S_0vT[i];
	ddS.S_0_nz++;
      }
    }
  }
  // check_m_vte(0);
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
	  if ( restart ) {
	    yap_quit("Inconsistency:  ddS.n_dt[%d][%d]=%d, ddS.c_dt[%d][%d]=%d\n",
		     i, t, (int)n, i, t, (int)ddS.c_dt[i][t]);
	  }
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
  if ( ddP.mu ) {
    ;   // do nothing
  } else if ( restart ) {
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
  if ( !ddP.mu ) {
    for (e=0; e<ddN.E; e++) 
      for (t=0; t<ddN.T; t++) 
	ddS.Cp_e[e] += ddS.cp_et[e][t];
  }
}

void check_cp_et() {
  int e, t;
  if ( ddP.mu )
    return;
  for (e=ddN.E-1; e>=0; e--) {
    for (t=0; t<ddN.T; t++) {
      int thisc = ddS.C_eDt[e][t];
      if ( e<ddN.E-1 )
	thisc += ddS.cp_et[e+1][t];
      if ( ((thisc>0) && (ddS.cp_et[e][t]==0)) 
	   || ((thisc>0) ^ (ddS.cp_et[e][t]>0)) 
	   || (ddS.cp_et[e][t]>thisc) ) {
	yap_quit("Inconsistency:  ddS.C_eDt[%d][%d]=%d, "
		 "ddS.cp_et[%d][%d]=%d\n",
		 e, t, (int)thisc, e, t, (int)ddS.cp_et[e][t]);
      }
    }
  }
}

void tca_reset_totals() {
  int e, i, t;
  /*
   *  reset not done for:
   *         ddS.n_dt, ddS.N_dT, ddS.c_dt, ddS.cp_et
   *         ddS.s_vte,  ddS.m_vte
   */
  memset((void*)ddS.M_Vte[0], 0, sizeof(ddS.M_Vte[t][e])*ddN.E*ddN.T);
  if ( !ddP.phi ) {
    memset((void*)ddS.S_Vte[0], 0, sizeof(ddS.S_Vte[t][i])*ddN.E*ddN.T);
    memset((void*)ddS.S_0vT, 0, sizeof(ddS.S_0vT[i])*ddN.W);
  }
  ddS.S_0 = 0;
  ddS.S_0_nz = 0;

  memset((void*)ddS.C_dT, 0, sizeof(ddS.C_dT[i])*ddN.D);
  memset((void*)ddS.C_e, 0, sizeof(ddS.C_e[i])*ddN.E);
  memset((void*)ddS.C_eDt[0], 0, sizeof(ddS.C_eDt[i][t])*ddN.E*ddN.T);
  if ( !ddP.mu ) 
    memset((void*)ddS.Cp_e, 0, sizeof(ddS.Cp_e[i])*ddN.E);

  /*
   *     beta part
   */
  for (e=ddN.E-1; e>=0; e--)
    for (i=0; i<ddN.W; i++) {
      for (t=0; t<ddN.T; t++) {
	ddS.M_Vte[t][e] += ddS.m_vte[i][t][e];
	if ( !ddP.phi ) ddS.S_Vte[t][e] += ddS.s_vte[i][t][e];
      } 
    }
  if ( !ddP.phi ) {
    for (i=0; i<ddN.W; i++) {
      for (t=0; t<ddN.T; t++) 
	ddS.S_0vT[i] += ddS.s_vte[i][t][0];
      assert(ddS.S_0vT[i]>=0);
      if ( ddS.S_0vT[i]>0 ) {
	ddS.S_0 += ddS.S_0vT[i];
	ddS.S_0_nz++;
      }
    }
  }
  /*
   *     alpha part
   */
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
  if ( !ddP.mu ) {
    for (e=0; e<ddN.E; e++) 
      for (t=0; t<ddN.T; t++) 
	ddS.Cp_e[e] += ddS.cp_et[e][t];
  }
}


void check_n_dt(int d) {
  int t;
  for (t=0; t<ddN.T; t++) {
    if ( ddS.n_dt[d][t]==0 && ddS.c_dt[d][t]>0 )
      assert(ddS.n_dt[d][t]>0 || ddS.c_dt[d][t]==0);
    if ( ddS.n_dt[d][t]>0 && ddS.c_dt[d][t]==0 )
      assert(ddS.n_dt[d][t]==0 || ddS.c_dt[d][t]>0);
  }
}

double global_sparsity(int k) {
  int v, e;
  double sparsity = 0;
  for (v=0; v<ddN.W; v++) {
    for (e=0; e<ddN.E; e++) {
      if ( ddS.m_vte[v][k][e]>0 ) {
        sparsity++;
        break;
      }
    }
  }
  return 1.0 - sparsity/ddN.W;
}

void check_m_vte(int e) {
  int v, t;
  if ( ddP.phi )
    return;
  for (v=0; v<ddN.W; v++) {
    for (t=0; t<ddN.T; t++) {
      int m = ddS.m_vte[v][t][e] +((e<ddN.E-1)?ddS.s_vte[v][t][e+1]:0);
      int s = ddS.s_vte[v][t][e];
      if ( m==0 && s>0 )
	assert(m>0 || ddS.s_vte[v][t][e]==0);
      if ( m>0 && s==0 )
	assert(m==0 || ddS.s_vte[v][t][e]>0);
    }
  }
  for (t=0; t<ddN.T; t++) {
    int total = 0;
    for (v=0; v<ddN.W; v++)
      total += ddS.m_vte[v][t][e];
    assert(total==ddS.M_Vte[t][e]);
    total = 0;
    for (v=0; v<ddN.W; v++)
      total += ddS.s_vte[v][t][e];
    assert(total==ddS.S_Vte[t][e]);
  }
  if ( e==0 ) {
    int tot=0;
    int totnz=0;
    for (v=0; v<ddN.W; v++) {
      int total = 0;
      for (t=0; t<ddN.T; t++) 
        total += ddS.s_vte[v][t][0];
      assert(total==ddS.S_0vT[v]);
      tot += ddS.S_0vT[v];
      if ( ddS.S_0vT[v]>0 )
        totnz++;
    }
    assert(tot==ddS.S_0);
    assert(totnz==ddS.S_0_nz);        
  }
}


