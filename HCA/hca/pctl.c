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
  for (p=ParA; p<=ParBeta; p++) {
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

  ddP.n_excludetopic = 0;
  ddP.excludetopic = NULL;
  ddP.bits_et = NULL;
  ddP.teststem = NULL;
  ddP.training = 0;
  ddP.memory = 0;
  ddP.phi = NULL;
  ddP.fixalpha = NULL;
  for (par=0; par<=ParBeta; par++) {
    ddT[par].samplerk = NULL;
    ddT[par].sampler = NULL;
    ddT[par].name = NULL;
    ddT[par].ptr = NULL;
    ddT[par].fix = 0;
    ddT[par].start = STARTCYCLES;
    ddT[par].cycles = 0;
    ddT[par].offset = 0;
  }
  ddT[ParAD].name = "ad";
  ddT[ParBDK].name = "bdk";
  ddT[ParA].name = "a";
  ddT[ParAW].name = "aw";
  ddT[ParA0].name = "a0";
  ddT[ParAW0].name = "aw0";
  ddT[ParB].name = "b";
  ddT[ParBW].name = "bw";
  ddT[ParB0].name = "b0";
  ddT[ParBW0].name = "bw0";
  ddT[ParAlpha].name = "alpha";
  ddT[ParBeta].name = "beta";
  ddT[ParA].ptr = &ddP.apar;
  ddT[ParAW].ptr = &ddP.awpar;
  ddT[ParA0].ptr = &ddP.a0;
  ddT[ParAW0].ptr = &ddP.aw0;
  ddT[ParB].ptr = &ddP.bpar;
  ddT[ParBW].ptr = &ddP.bwpar;
  ddT[ParB0].ptr = &ddP.b0;
  ddT[ParBW0].ptr = &ddP.bw0;
  ddT[ParAlpha].ptr = &ddP.alpha;
  ddT[ParBeta].ptr = &ddP.beta;
  ddT[ParAlpha].ptr = &ddP.alpha;
  ddT[ParAD].ptr = &ddP.ad;
  ddT[ParBDK].ptr = NULL;
  ddT[ParA].sampler = &sample_a;
  ddT[ParAW].sampler = &sample_aw;
  ddT[ParA0].sampler = &sample_a0;
  ddT[ParAW0].sampler =  &sample_aw0;
  ddT[ParB].sampler = &sample_b;
  ddT[ParBW].sampler = &sample_bw;
  ddT[ParB0].sampler = &sample_b0;
  ddT[ParBW0].sampler = &sample_bw0;
  ddT[ParAlpha].sampler = &sample_alpha;
  ddT[ParBeta].sampler = &sample_beta;
  ddT[ParAD].sampler = &sample_adk;
  ddT[ParBDK].samplerk = &sample_bdk;

  ddP.alpha = 1;
  ddP.betapr = NULL;
  ddP.betac = 1;
  ddP.beta = 0;
  ddP.PYalpha = H_HPDD;
  ddP.PYbeta = H_HPDD;
  ddP.apar = APAR;
  ddP.bpar = BPAR;
  ddP.awpar = AWPAR;
  ddP.bwpar = BWPAR;
  ddP.a0 = A0PAR;
  ddP.b0 = B0PAR;  
  ddP.aw0 = AW0PAR;
  ddP.bw0 = BW0PAR;  
  ddP.ad = APAR;
  ddP.bdk = NULL;  
  ddP.kbatch = 0;
  ddT[ParAlpha].cycles = DIRCYCLES;
  ddT[ParBeta].cycles = DIRCYCLES;
  ddT[ParB].cycles = BCYCLES;
  ddT[ParBDK].cycles = BCYCLES;
  ddT[ParB0].cycles = BCYCLES;
  ddT[ParBW].cycles = BCYCLES;
  ddT[ParBW0].cycles = BCYCLES;
  ddT[ParAD].cycles = ACYCLES;
  ddT[ParA].cycles = ACYCLES;
  ddT[ParA0].cycles = ACYCLES;
  ddT[ParAW].cycles = ACYCLES;
  ddT[ParAW0].cycles = ACYCLES;
  ddT[ParBeta].offset = 1;
  ddT[ParB0].offset = 1;
  ddT[ParBDK].offset = 0;
  ddT[ParBW0].offset = 1;
  ddT[ParA0].offset = 1;
  ddT[ParAW0].offset = 1;
  ddT[ParAD].offset = 2%BCYCLES;

  ddP.progiter = 5;
  ddP.progburn = 0;
  ddP.phiiter = 0;
  ddP.phiburn = 0;
  ddP.alphaiter = 0;
  ddP.alphaburn = 0;
  ddP.probiter = 0;
  ddP.probburn = 0;
  ddP.tprobiter = 0;
  ddP.tprobburn = 0;
  ddP.spiter = 0;
  ddP.spburn = 0;
  ddP.prditer = 15;
  ddP.prdburn = 5;
  ddP.lrsiter = 0;
  ddP.lrsburn = 0;
  ddP.mltiter = 15;
  ddP.mltburn = 5;
  ddP.cofile = NULL;
  ddP.queryiter = 0;
  ddP.query = NULL;
  ddP.qword = NULL;
  ddP.n_query = 0;

  ddP.Tinc = 5;
  ddP.Tcycle = 20;
  ddP.Tinit = 0;
  ddP.Tfree = -1;
  
  ddP.window = 0;
  ddP.window_cycle = 10;
  ddP.window_incr = 0;
  ddP.window_left = 0;
  ddP.window_right = 0;
 
  ddP.hold_every = 0;
  ddP.hold_dict = 0;
  ddP.hold_fraction = 0;
  ddP.docstats = NULL;
}

static char *mystem;
static char *mybuf;

static double readf(char *type) {
  char *par = readpar(mystem,type,mybuf,50);
  if ( par )
    return atof(par);
  return 0.0;
}
static int readi(char *type) {
  char *par = readpar(mystem,type,mybuf,50);
  if ( par )
    return atoi(par);
  return 0;
}
static double *readfv(char *type, int dim) {
  char *par = readpar(mystem,type,mybuf,dim*15+50);
  if ( par ) {
    double *vec = malloc(sizeof(vec[0])*dim);
    char *ptr;
    int t = 0;
    par += strspn(par, " ,");
    if ( !vec )
      yap_quit("Out of memory reading vector parameter '%s'\n", type);
    ptr = strtok(par," ,");
    while ( t<dim && ptr ) {
      vec[t] = atof(ptr);
      t++;
      ptr = strtok(NULL," ,");
    }
    if ( t<dim )
      yap_quit("Reading vector parameter '%s' only got %d/%d elements\n", 
	       type, t, dim);
    return vec;
  }
  return NULL;
}
static int *readiv(char *type, int dim) {
  char *par = readpar(mystem,type,mybuf,dim*15+50);
  if ( par ) {
    int *vec = malloc(sizeof(vec[0])*dim);
    char *ptr;
    int t = 0;
    par += strspn(par, " ,");
    if ( !vec )
      yap_quit("Out of memory reading vector parameter '%s'\n", type);
    ptr = strtok(par," ,");
    while ( t<dim && ptr ) {
      vec[t] = atoi(ptr);
      t++;
      ptr = strtok(NULL," ,");
    }
    if ( t<dim )
      yap_quit("Reading vector parameter '%s' only got %d/%d elements\n", 
	       type, t, dim);
    return vec;
  }
  return NULL;
}


void pctl_read(char *resstem, char *buf) {
  mystem = resstem;
  mybuf = buf;
  ddP.PYbeta = readi("PYbeta");
   if ( ddP.PYbeta ) {
     ddP.awpar = readf("aw");
     ddP.bwpar = readf("bw");
      if ( ddP.PYbeta!=H_PDP ) {
	ddP.aw0 = readf("aw0");
	ddP.bw0 = readf("bw0");
      }
   } else {
     ddP.beta = readf("beta");
   }
   ddP.PYalpha = readi("PYalpha");
   if ( ddP.PYalpha ) {
     ddP.apar = readf("a");
     ddP.bpar = readf("b");
     if ( ddP.PYalpha!=H_PDP ) {
	ddP.a0 = readf("a0");
	ddP.b0 = readf("b0");
      }
   } else {
     ddP.alpha = readf("alpha");
   }
   mybuf = malloc(ddN.T*15+100);
   if ( !mybuf ) 
     yap_quit("Out of memory reading 'bdk' in pctl_read()\n");
   ddP.bdk = readfv("bdk", ddN.T);
   if ( ddP.bdk!=NULL ) {
     ddP.ad = readf("ad");
   } else
     ddP.ad = 0;
   free(mybuf);
   mybuf = buf;
   ddP.n_excludetopic = readi("Nexcludetopic");
   if ( ddP.n_excludetopic>0 ) {
     int t, n_t;
     mybuf = malloc(ddP.n_excludetopic*10+100);
     ddP.excludetopic = readiv("excludetopic", ddP.n_excludetopic);
     /*  set the bit vector */
     n_t = ((ddN.T-1U)/32U+1U);
     ddP.bits_et = malloc(sizeof(ddP.bits_et[0])*n_t);
     for (t=0; t<n_t; t++) 
       ddP.bits_et[t] = 0;
     for (t=0; t<ddP.n_excludetopic; t++) {
       uint32_t x = ddP.excludetopic[t];
       ddP.bits_et[x/32U] |= (1U << (x%32U));
     }
     free(mybuf);
     mybuf = buf;
   }
}

/*
 *    fix pars based on data dims
 */
void pctl_dims() {
  if ( ddP.PYalpha==H_None ) {
    if ( ddP.alpha==0 ) {
      ddP.alpha = 0.05*ddN.NT/(ddN.DT*ddN.T);
    }
    if ( ddP.alpha< DIR_MIN )
      ddP.alpha = DIR_MIN*2.0;
    if ( ddP.alpha> DIR_MAX )
      ddP.alpha = DIR_MAX/2.0;
  }
  if ( ddP.PYbeta==H_None ) {
    if ( ddP.beta< DIR_MIN*ddN.W ) 
      ddP.beta = DIR_MIN*ddN.W*2.0;
    if ( ddP.beta>DIR_MAX*ddN.W ) 
      ddP.beta = DIR_MAX*ddN.W/2.0;
    if ( ddP.beta>DIR_TOTAL_MAX )
      ddP.beta = DIR_TOTAL_MAX/2.0;
  }
  if ( ddP.window>0 ) {
    if ( ddP.window>=ddN.DT )
      ddP.window = 0;
    ddP.window_left = 0;
    ddP.window_right = ddP.window;
  }
}

/*  needs to know ddN.T to work */
void pctl_fix(char *betafile, int ITER) {
  ddP.betapr = NULL;
  if ( ddP.ad==0 ) {
    ddT[ParAD].fix = 1;
  }
  if ( ddP.bdk!=NULL) {
    ddT[ParBDK].ptr = ddP.bdk;
    if ( ddT[ParBDK].fix==0 && ddP.kbatch==0 ) {
	if ( ddN.T>=20 )
      		ddP.kbatch = ddN.T/5;
	else
		ddP.kbatch = ddN.T/2;
    }
    if (  ddP.kbatch > ddN.T )
      ddP.kbatch = ddN.T;
  } else {
    ddT[ParBDK].fix = 1;   
  } 
  if ( ddP.bdk==NULL ) {
    ddT[ParAD].fix = 1;
  }
  if ( ddP.PYalpha==H_None ) {
    ddT[ParA].fix = 1;
    ddT[ParA0].fix = 1;
    ddT[ParB].fix = 1;
    ddT[ParB0].fix = 1;     
    ddP.alphaiter = 0;
    ddP.alphaburn = 0;
  } else {
    ddT[ParAlpha].fix = 1;
    if ( ddP.PYalpha==H_HDP )
      ddP.a0 = 0;
    if ( ddP.apar==0 )
      ddT[ParA].fix = 1;
    if ( ddP.a0==0 || ddP.PYalpha==H_PDP )
      ddT[ParA0].fix = 1;    
    if ( ddP.PYalpha==H_PDP )
      ddT[ParB0].fix = 1;    
  }
  if ( ddP.PYbeta==H_None ) {
    ddT[ParAW].fix = 1;
    ddT[ParAW0].fix = 1;
    ddT[ParBW].fix = 1;
    ddT[ParBW0].fix = 1;
  } else {
    if ( ddP.PYbeta==H_HDP ) 
      ddP.aw0 = 0;
    if ( ddP.awpar==0 )
      ddT[ParAW].fix = 1;
    if ( ddP.aw0==0 )
      ddT[ParAW0].fix = 1;    
    if ( ddP.PYbeta==H_PDP ) {
      ddT[ParAW0].fix = 1;    
      ddT[ParBW0].fix = 1;    
      ddP.beta = 1;
    }
    if ( ddP.PYbeta==H_HDP ) {
      ddP.beta = 1;
      if ( betafile==NULL )
	yap_quit("Option -Bhdp needs a '-u' too\n");
    }
    ddP.betac = 0;
    ddT[ParBeta].fix = 1;
    if ( betafile && ddP.PYbeta==H_HDP && ddP.aw0>0 ) {
      yap_quit("Only use '-u' when aw0==0\n");
    }
  }
  if ( ddP.phi!=NULL ) {
    /*
     *   PYtheta and beta are not used!!
     */
    ddT[ParAW].fix = 1;
    ddT[ParAW0].fix = 1;
    ddT[ParBW].fix = 1;
    ddT[ParBW0].fix = 1;
    ddT[ParBeta].fix = 1;
  }
  if ( ddP.fixalpha!=NULL ) {
    /*
     *   PYtheta and beta are not used!!
     */
    ddT[ParAlpha].fix = 1;
    ddT[ParA0].fix = 1;
    ddT[ParB0].fix = 1;
  }
  
  {
    enum ParType par;
    for (par=ParA; par<=ParBeta; par++) 
      ddT[par].offset =  ddT[par].offset %  ddT[par].cycles;
  }

  if ( ddP.query!=NULL && ddP.queryiter==0 )
    ddP.queryiter = 10;

  if ( ddP.lrsiter>0 ) {
    if ( ddP.lrsburn>=ddP.lrsiter )
      ddP.lrsiter = ddP.lrsburn+1;
    if ( ddP.lrsburn<1 )
      ddP.lrsburn = 1;
  } else
    ddP.lrsburn = 0;
  if ( ddP.mltiter>0 ) {
    if ( ddP.mltburn>=ddP.mltiter )
      ddP.mltiter = ddP.mltburn+1;
    if ( ddP.mltburn<1 )
      ddP.mltburn = 1;
  } else
    ddP.mltburn = 0;
  if ( ddP.prditer>0 ) {
    if ( ddP.prdburn>=ddP.prditer )
      ddP.prditer = ddP.prdburn+1;
    if ( ddP.prdburn<1 )
      ddP.prdburn = 1;
  } else
    ddP.prdburn = 0;

  if ( ddP.spiter==1 )
    ddP.spiter = 2;
  if ( ddP.probiter==1 )
    ddP.probiter = 2;
  if ( ddP.tprobiter==1 )
    ddP.tprobiter = 2;
  if ( ddP.phiiter==1 )
    ddP.phiiter = 2;
  if ( ddP.alphaiter==1 )
    ddP.alphaiter = 2;

  if ( ddP.Tinit==0 )
    ddP.Tinit = ddN.T;
  if ( ddP.Tfree<0 )
    ddP.Tfree = ITER;
}

int pctl_Tmax(int Tmax, int iter) {
  if ( Tmax<ddN.T && iter>(ddP.Tcycle*1.5) 
       && (ddP.Tcycle==1 || iter%ddP.Tcycle==0) ) {
    Tmax += ddP.Tinc;
    if ( Tmax>ddN.T )
      Tmax = ddN.T;
    if ( iter>=ddP.Tfree )
      Tmax = ddN.T;
  }
  return Tmax;
}
/*
 *   initialising or ddP.beta is changed
 */
void fixbeta(char *file, char *resstem) {
  int c;
  if ( ((ddP.PYbeta==H_PDP||ddP.PYbeta==H_HDP) && file) 
       || (ddP.PYbeta==H_None && !ddP.betapr) ) {
    /*
     *   only read on initialisation
     */
    ddP.betapr = dvec(ddN.W);
    if ( file ) {
      if ( strcmp(file,"file")==0 ) {
	for (c=0; c<ddN.W; c++)
	  ddP.betapr[c] = 0.5;
	for (c=0; c<ddN.N; c++) 
	  ddP.betapr[ddD.w[c]]++;
	for (c=0; c<ddN.W; c++)
	  ddP.betapr[c] /= ddN.N;
      } else if ( strcmp(file,"uniform")==0 ) {
	for (c=0; c<ddN.W; c++)
	  ddP.betapr[c] = 1.0/ddN.W;
      } else {
	double lastbeta = 0;
	assert(ddP.betac==0);
	read_dvec(file,ddN.W,ddP.betapr);
	//  normalise
	for (c=0; c<ddN.W; c++)
	  lastbeta += ddP.betapr[c];
	for (c=0; c<ddN.W; c++)
	  ddP.betapr[c] /= lastbeta;
      }
      if ( resstem ) {
	char *fname;
	fname = yap_makename(resstem,".beta");
	write_dvec(fname,ddN.W,ddP.betapr);
	free(fname);
      }
    }
  }
  if ( ddP.betac!=0 && ddP.PYbeta==H_None ) {
    assert(ddP.beta>0);
    ddP.betac = ddP.beta/ddN.W;
    for (c=0; c<ddN.W; c++)
      ddP.betapr[c] = ddP.betac;
  } else if ( ddP.betapr ) {
    double lastbeta = 0;
    for (c=0; c<ddN.W; c++)
      lastbeta += ddP.betapr[c];
    for (c=0; c<ddN.W; c++)
      ddP.betapr[c] *= ddP.beta/lastbeta;
  }
}


void pctl_report() {
  yap_message("PYbeta  = %d\n", (int)ddP.PYbeta);
  yap_message("beta  = %lf\n", ddP.beta);
  if ( ddP.betapr && ddP.betac==0 ) 
    yap_message("# beta from file\n");
  if ( ddP.PYbeta ) {
    yap_message("aw     = %lf\n", ddP.awpar);
    yap_message("bw     = %lf\n", ddP.bwpar);
    if ( ddP.PYbeta!=H_PDP ) {
      yap_message("aw0     = %lf\n", ddP.aw0);
      yap_message("bw0     = %lf\n", ddP.bw0);
    }
  }
  yap_message("PYalpha  = %d\n", (int)ddP.PYalpha);
  if ( ddP.PYalpha ) {
    yap_message("a     = %lf\n", ddP.apar);
    yap_message("b     = %lf\n", ddP.bpar);
    if ( ddP.PYalpha!=H_PDP ) {
      yap_message("a0     = %lf\n", ddP.a0);
      yap_message("b0     = %lf\n", ddP.b0);
    }
  } else {
    yap_message("alpha = %lf\n", ddP.alpha);
  }
  if ( ddP.bdk!=NULL ) {
    int t;
    yap_message("bdk =");
    for (t=0; t<ddN.T; t++) 
      yap_message(" %5lf", ddP.bdk[t]);
    yap_message("\n");
  }
  if ( ddP.bdk!=NULL )
    yap_message("ad  = %lf\n", ddP.ad);
  if ( ddP.n_excludetopic ) {
    int t;
    yap_message("excludetopic[%d] =", (int)ddP.n_excludetopic);
    for (t=0; t<ddP.n_excludetopic; t++) 
      yap_message(" %d", (int)ddP.excludetopic[t]);
    yap_message("\n");
  }
}

double pctl_gammaprior(double x) {
  return -x/PYP_CONC_PSCALE/ddN.W + (PYP_CONC_PSHAPE-1)*log(x);
}

/*
 *   generate parameter corresponding to index
 *   and return in *par and *k 
 *        note bdk has K values bdk[k]
 *        all other pars are single valued
 *   return 1 if found OK, else return 0 if no more
 */
int pctl_par_iter(int index, int iter, enum ParType *par, int *k) {
  enum ParType p;
  for (p=ParA; p<=ParBeta; p++) {
    if (  !ddT[p].fix && ddT[p].ptr
          && iter>ddT[p].start
	  && iter%ddT[p].cycles==ddT[p].offset ) {
      if ( p==ParBDK ) {
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
  int p, index;
  struct pst_data pd;
#ifdef H_THREADS
  pthread_t thread[procs];
#endif
  
  /*
   *  first, create docstats if needed
   */
  ddP.docstats = NULL;
  for (index=0; index<100000; index++) {
    int k;
    enum ParType par;
    if ( pctl_par_iter(index, iter, &par, &k) && par==ParBDK ) {
      ddP.docstats = dmi_bstore(&ddM);
      break;
    }
  }
  index = 0;
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
  if ( ddP.docstats ) {
    dmi_freebstore(&ddM,ddP.docstats);
    ddP.docstats = NULL;
  }
}

/*
 *    i = word index in z[], w[] etc.
 *    if its in training set, cannot be hold out
 *    otherwise compute hold out
 */
int pctl_hold(int i) {
  if ( i>=ddN.NT ) { 
    int d = ddD.d[i];
    int starti = ddD.NdTcum[d];
    if ( ddP.hold_dict ) {
      if ( ((ddD.w[i]+1)%ddP.hold_dict)==0 )
	return 1;
    } else if ( ddP.hold_every ) {
      if ( (i-starti+1)%ddP.hold_every==0 )
	return 1;
    } else {
      if ( i-starti>ddD.NdT[d]*ddP.hold_fraction )
	return 1;
    }
  }
  return 0;
}

void pctl_samplereport() {
  enum ParType par;
  yap_message("    sampling pars:");
  for (par=ParA; par<=ParBeta; par++) {
    if (  !ddT[par].fix )
      yap_message(" %s(%d),", ddT[par].name, ddT[par].cycles);
  }
  yap_message("\n");
}

void pctl_update(int iter) {
  enum ParType par;
  int start=1;
  yap_message("Pars:  ");
  for (par=ParA; par<=ParBeta; par++) {
    if (  !ddT[par].fix && iter>ddT[par].start ) {
      if ( !start ) 
	yap_message(", ");
      yap_message("%s=%lf", ddT[par].name, *ddT[par].ptr);
      start = 0;
    }
  }  
  yap_message("\n");
}

static void printpar(FILE *fp, enum ParType par) {
  if ( !ddT[par].fix ) 
    fprintf(fp, "#  %s was sampled every %d major cycles\n", 
	    ddT[par].name, ddT[par].cycles);
  fprintf(fp, "%s = %lf\n", ddT[par].name, *ddT[par].ptr);
}

void pctl_print(FILE *fp) {
  printpar(fp,ParBeta);
  fprintf(fp, "PYbeta  = %d\n", (int)ddP.PYbeta);
  if ( ddP.PYbeta ) {
    printpar(fp,ParAW); printpar(fp,ParBW);
    if ( ddP.PYbeta !=H_PDP ) {
      printpar(fp,ParAW0); printpar(fp,ParBW0);
    }
  }
  fprintf(fp, "PYalpha  = %d\n", (int)ddP.PYalpha);
  if ( ddP.PYalpha ) {
    printpar(fp,ParA); printpar(fp,ParB);
    if ( ddP.PYalpha!=H_PDP ) {
      printpar(fp,ParA0); printpar(fp,ParB0);
    }
  } else {
    printpar(fp,ParAlpha);
  }
  if ( ddP.bdk!=NULL ) {
    int t;
    if ( !ddT[ParBDK].fix ) 
      fprintf(fp, "#  %s was sampled every %d major cycles in batches of %d\n", 
	      ddT[ParBDK].name, ddT[ParBDK].cycles, ddP.kbatch);
    fprintf(fp, "bdk =");
    for (t=0; t<ddN.T; t++) 
      fprintf(fp, " %5lf", ddP.bdk[t]);
    fprintf(fp, "\n");
    fprintf(fp, "# %s S-table had bounds N=%d and T=%d\n",
	    ddT[ParAD].name, ddC.SD->maxN, ddC.SD->maxM);
    printpar(fp,ParAD);
  }
  if ( ddP.n_excludetopic ) {
    int t;
    fprintf(fp, "Nexcludetopic = %d\n", ddP.n_excludetopic);
    fprintf(fp, "excludetopic =");
    for (t=0; t<ddP.n_excludetopic; t++) 
      fprintf(fp, " %d", ddP.excludetopic[t]);
    fprintf(fp, "\n");
  }
}

void pctl_free() {
  if ( ddP.bdk!=NULL )
    free(ddP.bdk);
  if ( ddP.phi ) {
    free(ddP.phi[0]); free(ddP.phi);
  }
  if ( ddP.fixalpha ) {
    free(ddP.fixalpha);
  }
  if ( ddP.n_excludetopic ) {
    free(ddP.excludetopic);
    free(ddP.bits_et);
  }
  if ( ddP.qword ) {
    free(ddP.qword);
    free(ddP.qposn);
    free(ddP.query);
    free(ddP.qid);
  }
  if ( ddP.betapr )
    free(ddP.betapr);
  
}
