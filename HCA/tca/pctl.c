/*
 * Hyperparameters and their control
 * Copyright (C) 2013-2014 Wray Buntine
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

#include "pctl.h"
#include "util.h"
#include "yap.h"
#include "sample.h"
#include "data.h"
#ifdef H_THREADS
#include <pthread.h>
#endif
#include "atomic.h"

enum ParType findpar(char *name) {
  enum ParType p;
  for (p=ParAM; p<=ParBB; p++) {
    if ( strcmp(name,ddT[p].name)==0 )
      return p;
  }
  return ParNone;
}

int pctl_training(int D) {
  int train = 0;
  if ( ddP.training==0 )
    train = D-ddN.TEST;
  else if ( ddP.training<=D-ddN.TEST )
    train = ddP.training;
  else {
    train = D-ddN.TEST;
  }
  return train;
}

void pctl_init() {
  enum ParType par;

  for (par=0; par<=ParBB; par++) {
    ddT[par].sampler = NULL;
    ddT[par].samplerk = NULL;
    ddT[par].name = NULL;
    ddT[par].ptr = NULL;
    ddT[par].fix = 0;           
    ddT[par].start = 10;
    ddT[par].cycles = 0;
    ddT[par].offset = 0;
  }
  ddT[ParAM].name = "am";
  ddT[ParAT].name = "at";
  ddT[ParAB].name = "ab";
  ddT[ParAP0].name = "ap0";
  ddT[ParAP1].name = "ap1";
  ddT[ParBM0].name = "bm0";
  ddT[ParBM1].name = "bm1";
  ddT[ParBT].name = "bt";
  ddT[ParBB].name = "bb";
  ddT[ParB0P].name = "b0p";
  ddT[ParB0M].name = "b0m";
  ddT[ParBP0].name = "bp0";
  ddT[ParBP1].name = "bp1";

  ddT[ParAT].ptr = &ddP.a_theta;
  ddT[ParAM].ptr = &ddP.a_mu;
  ddT[ParAB].ptr = &ddP.a_burst;
  ddT[ParAP0].ptr = &ddP.a_phi0;
  ddT[ParAP1].ptr = &ddP.a_phi1;
  ddT[ParB0P].ptr = &ddP.b_phi0;
  ddT[ParB0M].ptr = &ddP.b_mu0;
  ddT[ParBT].ptr = &ddP.b_theta;
  ddT[ParBB].ptr = &ddP.b_burst;
  ddT[ParBM0].ptr = NULL;
  ddT[ParBM1].ptr = NULL;
  ddT[ParBP0].ptr = NULL;
  ddT[ParBP1].ptr = NULL;
  
  ddT[ParAT].sampler = &sample_at;
  ddT[ParAM].sampler = &sample_am;
  ddT[ParAB].sampler = &sample_ab;
  ddT[ParAP0].sampler = &sample_ap0;
  ddT[ParAP1].sampler = &sample_ap1;
  ddT[ParB0P].sampler = NULL;
  ddT[ParB0M].sampler = NULL;
  ddT[ParBM0].sampler = &sample_bm0;
  ddT[ParBM1].sampler = &sample_bm1;
  ddT[ParBT].sampler = &sample_bt;
  ddT[ParBB].sampler = &sample_bb;
  ddT[ParBP0].sampler = &sample_bp0;
  ddT[ParBP1].samplerk = &sample_bp1;

  ddT[ParAB].fix = 1;
  ddT[ParB0P].fix = 1;
  ddT[ParB0M].fix = 1;
  // ddT[ParBP0].fix = 1;
  //WRAY  sampling BM0 causes weird memory bugs
  //  ddT[ParBM0].fix = 1;
  ddT[ParAP0].fix = 1;
  ddT[ParAP1].fix = 1;

  ddP.a_mu = APAR;
  ddP.a_theta = APAR;
  ddP.a_phi0 = APAR;
  ddP.a_phi1 = AWPAR;
  ddP.a_burst = AWPAR;
  ddP.b_mu = NULL;
  ddP.b_mu0 = B0PAR;
  ddP.b_theta = BPAR;
  ddP.b_burst = 0;
  ddP.b_phi = NULL;
  ddP.b_phi0 = B0PAR;

  ddT[ParBM0].cycles = BCYCLES;
  ddT[ParBM1].cycles = BCYCLES;
  ddT[ParBP0].cycles = BCYCLES;
  ddT[ParBP1].cycles = BCYCLES;
  ddT[ParBT].cycles = BCYCLES;
  ddT[ParBB].cycles = BCYCLES;
  ddT[ParBB].offset = 1%BCYCLES;
  ddT[ParBM0].offset = 1%BCYCLES;
  ddT[ParBM1].offset = 1%BCYCLES;
  
  ddT[ParAT].cycles = ACYCLES;
  ddT[ParAB].cycles = ACYCLES;
  ddT[ParAM].cycles = ACYCLES;
  ddT[ParAP0].cycles = ACYCLES;
  ddT[ParAP1].cycles = ACYCLES;
  ddT[ParAB].offset = 1%ACYCLES;
  ddT[ParAM].offset = 2%ACYCLES;
  ddT[ParAP0].offset = 3%ACYCLES;
  ddT[ParAP1].offset = 3%ACYCLES;

  ddP.phiiter = 0;
  ddP.phiburn = 0;
  ddP.muiter = 0;
  ddP.muburn = 0;
  ddP.progiter = 2;
  ddP.progburn = 1;
  ddP.mltiter = 15;
  ddP.mltburn = 5;
  ddP.cofile = NULL;
  ddP.teststem = NULL;
  ddP.training = 0;
  ddP.back = UINT32_MAX;
  ddP.kbatch = 0;
  
  ddP.maxN = 20000;
  ddP.maxM = 1000;

  ddP.hold_every = 0;
  ddP.hold_dict = 0;
  ddP.hold_fraction = 0;

  ddP.mu = NULL;
  ddP.phi = NULL;
}

static char *mystem;
static char *mybuf;

static double readf(char *type) {
  char *par = readpar(mystem,type,mybuf,50);
  if ( par )
    return atof(par);
  return 0.0;
}
static double *readfv(char *type, int dim) {
  char *par = readpar(mystem,type,mybuf,dim*15+strlen(type)+3);
  if ( par ) {
    double *vec = malloc(sizeof(vec[0])*dim);
    char *ptr;
    int t = 0;
    par += strspn(par, " ,");
    if ( !vec )
      yap_quit("Out of memory reading vector parameter '%s'\n", type);
    while ( t<dim && (ptr=strsep(&par," ,")) ) {
      vec[t] = atof(ptr);
      t++;
    }
    if ( t<dim )
      yap_quit("Reading vector parameter '%s' only got %d/%d elements\n", 
	       type, t, dim);
    return vec;
  }
  return NULL;
}


void pctl_read(char *resstem, char *buf) {
  int i;
  mybuf = malloc(ddN.T*15+10);
  if ( !mybuf )
    yap_quit("Out of memory in pctl_read()\n");
  mystem = resstem;
  ddP.a_mu = readf("am");
  ddP.a_theta = readf("at");
  ddP.a_phi0 = readf("ap0");
  ddP.a_phi1 = readf("ap1");
  ddP.a_burst = readf("ab");
  ddP.b_mu = (double*)malloc(sizeof(*ddP.b_mu)*ddN.E);
  ddP.b_mu[0] = readf("bm0");  
  if ( ddN.E>1 ) {
    ddP.b_mu[1] = readf("bm1");
    for (i=2; i<ddN.E; i++)
      ddP.b_mu[i] = ddP.b_mu[1];
  }
  ddP.b_mu0 = readf("b0m");
  ddP.b_theta = readf("bt");
  ddP.b_burst = readf("bb");
  ddP.b_phi0 = readf("b0p");
  ddP.b_phi = (double**)malloc(sizeof(*ddP.b_phi)*ddN.E);
  ddP.b_phi[0] = (double*)malloc(sizeof(*ddP.b_phi[0])*ddN.T);
  if ( !ddP.b_phi || !ddP.b_phi[0] )
    yap_quit("Out of memory in pctl_read()\n");
  ddP.b_phi[0][0] = readf("bp0");
  for (i=1; i<ddN.T; i++)
      ddP.b_phi[0][i] = ddP.b_phi[0][0];
  if ( ddN.E>1 ) {
    ddP.b_phi[1] = readfv("bp1", ddN.T);
    for (i=2; i<ddN.E; i++)
      ddP.b_phi[i] = ddP.b_phi[1];
  }
  free(mybuf);
}

void pctl_loadphi(char *resstem) {
  int e;
  char ebuf[10];
  char *fname;
  ddP.phi = malloc(sizeof(ddP.phi[0])*ddN.E);
  if ( ddN.E>1000 )
    yap_quit("Compiled topic size bound too small in pctl_loadphi()\n");
  for (e=0; e<ddN.E; e++) {
    ddP.phi[e] = fmat(ddN.W, ddN.T);
    if (  !ddP.phi[e] ) 
      yap_quit("Cannot allocate memory in pctl_loadphi()\n");
    sprintf(ebuf,".phi%03d",e);
    fname = yap_makename(resstem, ebuf);
    read_fmat(fname, ddN.E, ddN.T, ddP.mu);
    free(fname);
  }
}

void pctl_loadmu(char *resstem) {
  char *fname = yap_makename(resstem, ".mu");
  ddP.mu = fmat(ddN.E, ddN.T);
  if ( !ddP.mu ) 
    yap_quit("Cannot allocate memory in pctl_loadmu()\n");
  read_fmat(fname, ddN.E, ddN.T, ddP.mu);
  free(fname);
}


void pctl_free() {
  free(ddP.b_mu);
  ddP.b_mu = NULL;
  free(ddP.b_phi[0]);
  if ( ddN.E>1 ) 
    free(ddP.b_phi[1]);
  free(ddP.b_phi);
  ddP.b_phi = NULL;
  if ( ddP.mu ) {
    free(ddP.mu[0]); free(ddP.mu);
  }
  if ( ddP.phi ) {
    int e;
    for (e=0; e<ddN.E; e++) {
      free(ddP.phi[e][0]); free(ddP.phi[e]);
    } 
    free(ddP.phi);
  }
}

void pctl_fix(int ITER) {
  if ( ddP.a_mu==0 )
    ddT[ParAM].fix = 1;
  if ( ddP.a_theta==0 )
    ddT[ParAT].fix = 1;
  if ( ddP.a_phi0==0 )
    ddT[ParAP0].fix = 1;
  if ( ddP.a_phi1==0 )
    ddT[ParAP1].fix = 1;
  if ( ddP.b_burst==0 ) {
    ddT[ParAB].fix = 1;
    ddT[ParBB].fix = 1;
  }
  if ( ddP.b_phi==NULL ) {
    int i;
    ddP.b_phi = (double**)malloc(sizeof(*ddP.b_phi)*ddN.E);
    if ( !ddP.b_phi )
      yap_quit("Out of memory in pctl_fix()\n");
    ddP.b_phi[0] = dvec(ddN.T);
    for (i=0; i<ddN.T; i++) {
      ddP.b_phi[0][i] = BWPAR;
    }
    if ( ddN.E>1 ) {
      ddP.b_phi[1] = dvec(ddN.T);
      for (i=0; i<ddN.T; i++) {
        ddP.b_phi[1][i] = BWPAR;
      }
      for (i=2; i<ddN.E; i++)
        ddP.b_phi[i] = ddP.b_phi[1];
    }
  }
  if ( ddP.back==UINT32_MAX )
    ddP.back = ddN.E+1;
  if ( ddN.E==1 ) {
    ddT[ParBP1].fix = 1;
    ddT[ParBM1].fix = 1;
  }
  if ( ddP.b_mu==NULL ) {
    int i;
    ddP.b_mu = dvec(ddN.E);
    for (i=0; i<ddN.E; i++) {
      ddP.b_mu[i] = BPAR;
    }
  }
  ddT[ParBM0].ptr = &ddP.b_mu[0];
  ddT[ParBP0].ptr = ddP.b_phi[0];
  if ( ddN.E>1 ) {
    ddT[ParBM1].ptr = &ddP.b_mu[1];
    ddT[ParBP1].ptr = ddP.b_phi[1];
  }

  if ( ddP.mu ) {
    ddP.muiter = 0;
    ddT[ParAM].fix = 1; 
    ddT[ParBM1].fix = 1; 
    ddT[ParBM0].fix = 1; 
    ddT[ParB0M].fix = 1;
  }
  if ( ddP.phi ) {
    ddP.phiiter = 0;
    ddT[ParAP0].fix = 1; 
    ddT[ParAP1].fix = 1; 
    ddT[ParBP0].fix = 1; 
    ddT[ParBP1].fix = 1; 
    ddT[ParB0P].fix = 1; 
  }

  {
    enum ParType par;
    for (par=ParAM; par<=ParBB; par++) 
      if ( ddT[par].fix==0 )
	ddT[par].offset =  ddT[par].offset %  ddT[par].cycles;
  }

  if ( ddP.mltiter>0 ) {
    if ( ddP.mltburn>=ddP.mltiter )
      ddP.mltiter = ddP.mltburn+1;
    if ( ddP.mltburn<1 )
      ddP.mltburn = 1;
  } else
    ddP.mltburn = 0;

  if ( ddP.maxM>ddP.maxN )
    ddP.maxM = ddP.maxN;

  if ( ddT[ParBP1].fix==0 ) {
    /*
     *  adjust kbatch
     */
    if ( ddP.kbatch==0 ) {
      if ( ddN.T>=20 )
        ddP.kbatch = ddN.T/5;
      else
        ddP.kbatch = ddN.T/2;
    }
    if (  ddP.kbatch > ddN.T )
      ddP.kbatch = ddN.T;
  }
  if ( ddP.phiiter==1 )
    ddP.phiiter = 2;
  if ( ddP.muiter==1 )
    ddP.muiter = 2;
}

void pctl_report() {
  int i;
  yap_message("am = %lf\n", ddP.a_mu);
  yap_message("bm0 = %lf\n", ddP.b_mu[0]);
  if ( ddN.E>1 ) 
    yap_message("bm1 = %lf\n", ddP.b_mu[1]);
  yap_message("b0m = %lf\n", ddP.b_mu0);
  yap_message("at = %lf\n", ddP.a_theta);
  yap_message("bt = %lf\n", ddP.b_theta);
  yap_message("b0p = %lf\n", ddP.b_phi0);
  yap_message("ap0 = %lf\n", ddP.a_phi0);
  yap_message("ap1 = %lf\n", ddP.a_phi1);
  yap_message("bp0 = %lf", ddP.b_phi[0][0]);
  if ( ddN.E>1 ) {
    yap_message("\nbp1 = ");
    for (i=0; i<ddN.T; i++)
      yap_message("%lf,", ddP.b_phi[1][i]);
  }
  yap_message("\n");
  yap_message("ab = %lf\n", ddP.a_burst);
  yap_message("bb = %lf\n", ddP.b_burst);
}

double pctl_gammaprior(double x) {
  return -x/PYP_CONC_PSCALE/ddN.W + (PYP_CONC_PSHAPE-1)*log(x);
}

/*
 *   generate parameter corresponding to index
 *   and return in *par and *k 
 *        note b_phi has K values b_phi[1][k]
 *        all other pars are single valued
 *   return 1 if found OK, else return 0 if no more
 */
int pctl_par_iter(int index, int iter, enum ParType *par, int *k) {
  enum ParType p;
  for (p=ParAM; p<=ParBB; p++) {
    if (  !ddT[p].fix && ddT[p].ptr
          && iter>ddT[p].start
	  && iter%ddT[p].cycles==ddT[p].offset ) {
      if ( p==ParBP1 ) {
        if ( index<ddP.kbatch ) {
          *par = p;
          *k = (iter*ddP.kbatch+index)%ddN.T;
          return 1;
        }
        index -= ddP.kbatch;
      } else {
        if ( index==0 ) {
          *par = p;
          *k = -1;
          return 1;
        }
        index--;
      }
    }
  }
  return 0;
}
struct pst_data {
  int iter;
  int *index;  /*  shared location to get index */
};
static void *pctl_sample_thread(void *pin) {
  struct pst_data *pd=(struct pst_data *)pin;
  double startlike = 0;
  int k, index;
  enum ParType par;
  while ( 1 ) {
    index = atomic_incr(*pd->index) - 1;
    if ( pctl_par_iter(index, pd->iter, &par, &k) ) {
      if ( verbose>2 ) {
	/*  fetching likelihood very expensive!! */
        startlike = likelihood();
        if ( k<0 )
          yap_message("sample_%s", ddT[par].name);
        else
          yap_message("sample_%s[%d]", ddT[par].name, k);
        yap_message(" (pre): %s=%lf, lp=%lf\n",
                    ddT[par].name, ddT[par].ptr[k<0?0:k], startlike);
      }
      if ( k<0 )
        (*ddT[par].sampler)(ddT[par].ptr);
      else
        (*ddT[par].samplerk)(ddT[par].ptr,k);
      if ( verbose>2 ) {
        double endlike = likelihood();
        if ( k<0 )
          yap_message("sample_%s", ddT[par].name);
        else
          yap_message("sample_%s[%d]", ddT[par].name, k);
        yap_message(" (pre): %s=%lf, lp=%lf\n",
                    ddT[par].name, ddT[par].ptr[k<0?0:k], endlike);
        if ( pd->iter>50 && (endlike-startlike)/ddN.NT>1 ) {
          yap_quit("Sampler failed iter=%d due to huge decrease of %lf!\n",
                   pd->iter, (endlike-startlike)/ddN.NT);
        }
      }
    } else
      break;
  }
  return NULL;
}

void pctl_sample(int iter, int procs) {
#ifdef H_THREADS
  int p;
#endif
  int index = 0;
  struct pst_data pd;
#ifdef H_THREADS
  pthread_t thread[procs];
#endif
  pd.index = &index;
  pd.iter = iter;
#ifdef H_THREADS
  if ( procs>1 ) {
      for (p = 0 ; p < procs ; p++){ 
        if ( pthread_create(&thread[p],NULL,pctl_sample_thread,(void*)&pd) != 0) {
          yap_message("pctl_sample() thread failed %d\n",p+1 );
        }
      }
      //waiting for threads to finish
      for (p = 0; p < procs; p++){
        pthread_join(thread[p], NULL);
      }
  } else
#endif
  pctl_sample_thread((void*)&pd);
}

/*
 *    i = word index in z[], w[] etc.
 *    if its in training set, cannot be hold out
 *    otherwise compute hold out
 */
int pctl_hold(int i) {
  if ( i>=ddN.NT ) { 
    int d = ddD.d[i];
    int starti = ddD.N_dTcum[d];
    if ( ddP.hold_dict ) {
      if ( ((ddD.w[i]+1)%ddP.hold_dict)==0 )
	return 1;
    } else if ( ddP.hold_every ) {
      if ( (i-starti+1)%ddP.hold_every==0 )
	return 1;
    } else {
      if ( i-starti>ddD.N_dT[d]*ddP.hold_fraction )
	return 1;
    }
  }
  return 0;
}

void pctl_samplereport() {
  enum ParType par;
  yap_message("    sampling pars:");
  for (par=ParAM; par<=ParBB; par++) {
    if (  !ddT[par].fix ) {
      if ( par==ParBP1 ) 
        yap_message(" %sX%d(%d),", ddT[par].name, ddP.kbatch,ddT[par].cycles);
      else 
        yap_message(" %s(%d),", ddT[par].name, ddT[par].cycles);
    }
  }
  yap_message("\n");
}

void pctl_update(int iter) {
  enum ParType par;
  for (par=ParAM; par<=ParBB; par++) {
    if (  !ddT[par].fix && iter>ddT[par].start ) {
      yap_message(", %s=%lf", ddT[par].name, *ddT[par].ptr);
    }
  }  
  yap_message("\n");
}

static void printpar(FILE *fp, enum ParType par) {
  if ( ddT[par].ptr==NULL )
    return;
  if ( !ddT[par].fix ) 
    fprintf(fp, "#  %s was sampled every %d major cycles\n", 
	    ddT[par].name, ddT[par].cycles);
  fprintf(fp, "%s = %lf\n", ddT[par].name, *ddT[par].ptr);
}

void pctl_print(FILE *fp) {
  int i;
  printpar(fp,ParAM);
  printpar(fp,ParBM0);
  printpar(fp,ParBM1);
  printpar(fp,ParB0M);
  printpar(fp,ParAT);
  printpar(fp,ParBT);
  printpar(fp,ParAB);
  printpar(fp,ParBB);
  printpar(fp,ParAP0);
  printpar(fp,ParAP1);
  printpar(fp,ParB0P);
  if ( !ddT[ParBP0].fix ) 
    fprintf(fp, "#  %s was sampled every %d major cycles\n", 
	    ddT[ParBP0].name, ddT[ParBP0].cycles);
  fprintf(fp, "bp0 = %lf\n", ddP.b_phi[0][0]);
  if ( !ddT[ParBP1].fix ) 
    fprintf(fp, "#  %s was sampled every %d major cycles\n", 
	    ddT[ParBP1].name, ddT[ParBP1].cycles);
  if ( ddN.E>1 ) {
    fprintf(fp, "bp1 = ");
    for (i=0; i<ddN.T; i++)
      fprintf(fp, "%lf,", ddP.b_phi[1][i]);
    fprintf(fp, "\n");
  }
}
