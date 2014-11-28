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
  ddP.theta = NULL;
  ddP.phi = NULL;
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
  ddT[ParNGAlpha].name = "NGalpha"; 
  ddT[ParNGBeta].name = "NGbeta";  
  ddT[ParAlpha].name = "alphatot"; 
  ddT[ParBeta].name = "betatot";  
  ddT[ParA].ptr = &ddP.apar;
  ddT[ParAW].ptr = &ddP.awpar;
  ddT[ParA0].ptr = &ddP.a0;
  ddT[ParAW0].ptr = &ddP.aw0;
  ddT[ParB].ptr = &ddP.bpar;
  ddT[ParBW].ptr = &ddP.bwpar;
  ddT[ParB0].ptr = &ddP.b0;
  ddT[ParBW0].ptr = &ddP.bw0;
  ddT[ParNGAlpha].ptr = ddP.NGalpha;
  ddT[ParNGBeta].ptr = ddP.NGbeta;
  ddT[ParAlpha].ptr = &ddP.alphatot;
  ddT[ParBeta].ptr = &ddP.betatot;
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
  ddT[ParNGAlpha].samplerk = &sample_NGalpha;
  ddT[ParNGBeta].samplerk = &sample_NGbeta;
  ddT[ParAD].sampler = &sample_adk;
  ddT[ParBDK].samplerk = &sample_bdk;

  ddP.alphatot = 0;
  ddP.alphac = 0;
  ddP.NGalpha = NULL;
  ddP.NGbeta = NULL;
  ddP.alphapr = NULL;
  ddP.betapr = NULL;
  ddP.betac = 0;
  ddP.betatot = 0;
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
  ddT[ParNGAlpha].cycles = DIRCYCLES;
  ddT[ParNGBeta].cycles = DIRCYCLES;
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
  ddT[ParNGBeta].offset = 1;
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
  ddP.mergeiter = 0;
  ddP.mergeinit = 0;
  ddP.mergemin = 0;
  ddP.mergebest = 1;

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
  ddP.hold_all = 0;
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
  char *par = readpar(mystem,type,mybuf,dim*20+50);
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
  char *par = readpar(mystem,type,mybuf,dim*20+50);
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

/*
 *  reads parameters, not dimensions
 */
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
    if ( ddP.PYbeta!=H_HPDD )
      ddP.betatot = 1.0;
    else
      ddP.betatot = 0;
    ddP.betac = 0.0;
  } else {
    ddP.betatot = readf("betatot");
  }
  ddP.PYalpha = readi("PYalpha");
  if ( ddP.PYalpha ) {
    /*
     *  its a Pitman-Yor
     */
    ddP.apar = readf("a");
    ddP.bpar = readf("b");
    if ( ddP.PYalpha!=H_PDP ) {
      ddP.a0 = readf("a0");
      ddP.b0 = readf("b0");
    }
    if ( ddP.PYalpha!=H_HPDD && ddP.PYalpha!=H_NG )
      ddP.alphatot = 1.0;
    else
      ddP.alphatot = 0;
    ddP.alphac = 0.0;
  } else {
    ddP.NGalpha = readfv("NGalpha", ddN.T);
    if ( ddP.NGalpha ) {
      int t;
      double tot = 0;
      /*
       *  its a normalised gamma
       */
      ddP.NGbeta = readfv("NGbeta", ddN.T);
      if ( !ddP.NGbeta ) 
        yap_quit("Cannot read 'NGbeta' when 'NGalpha' exists in '%s.par'\n",
		 resstem);
      for (t=0; t<ddN.T; t++) tot += ddP.NGbeta[t];
      tot /= ddN.T;
      for (t=0; t<ddN.T; t++) ddP.NGbeta[t] /= tot;
      ddP.NGbetamin = dmin(ddN.T, ddP.NGbeta);
      ddP.NGbetamax = dmax(ddN.T, ddP.NGbeta);
    } else {
      /*
       *  its a Dirichlet
       */
      ddP.alphatot = readf("alphatot");
    }
  }
  mybuf = malloc(ddN.T*20+100);
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
 *   default alpha values for LDA
 */
static double pctl_alphacinit() {
  return 0.05*ddN.NT/(ddN.DT*ddN.T);
}
static double pctl_alpharange(double alphac) {
  if ( alphac< DIR_MIN ) 
    alphac = DIR_MIN;
  if ( alphac>DIR_MAX ) 
    alphac = DIR_MAX;
  if ( alphac>DIR_TOTAL_MAX/ddN.T )
    alphac = DIR_TOTAL_MAX/ddN.T;
  return alphac;
}

/*
 *    fix pars based on data dims ..
 *    called with pctl_fixbeta() at initialisation
 */
void pctl_dims() {
  if ( ddP.bdk!=NULL) {
    ddT[ParBDK].ptr = ddP.bdk;
  }
 if ( ddT[ParBDK].fix==0 || ddT[ParNGAlpha].fix==0 || ddT[ParNGBeta].fix==0 ) { 
    /*
     *    set kbatch for sampling
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
  if ( ddP.PYalpha==H_NG ) {
    /*  
     * initialisation sets up values like LDA default
     */
    int t;
    if ( !ddP.NGbeta ) {
      ddP.NGbeta = malloc(ddN.T*sizeof(*ddP.NGbeta));
      for (t=0; t<ddN.T; t++)
        ddP.NGbeta[t] = 1.0/ddN.T;
      ddP.NGbetamin = 1.0/ddN.T;
      ddP.NGbetamax = 1.0/ddN.T;
    }
    if ( !ddP.NGalpha ) {
      double alphac = pctl_alpharange(pctl_alphacinit());
      ddP.NGalpha = malloc(ddN.T*sizeof(*ddP.NGalpha));
      for (t=0; t<ddN.T; t++)
        ddP.NGalpha[t] = alphac;
    }
    ddP.PYalpha = H_None;
  }
  if ( ddP.PYalpha==H_None && !ddP.NGalpha ) {
    /*
     *  reset .alphatot according to constraints;
     *  if.alphatot== then .alphac==0.05*(NT/DT)/T;
     */
    double alphain = ddP.alphatot/ddN.T;
    double alphac = alphain;
    if ( alphac==0 )
      ddP.alphac = pctl_alphacinit();
    ddP.alphac = pctl_alpharange(ddP.alphac);
    if ( verbose>=1 && alphain!=alphac ) {
      yap_message("alpha changed from %lf to %lf due to Dirichlet constrains\n",
		  alphain, alphac);
      ddP.alphatot = alphac*ddN.T;
    }
  } 
  if ( ddP.PYbeta==H_None ) {
    /*
     *  reset .betatot according to constraints;
     *  if.betatot== then .betac==DIR_MIN;
     */
    double betain = ddP.betatot/ddN.W;
    double betac = betain;
     if ( betac==0 ) 
      betac = DIR_MIN*10;
     if ( betac< DIR_MIN ) 
      betac = DIR_MIN;
    if ( betac>DIR_MAX ) 
      betac = DIR_MAX;
    if ( betac>DIR_TOTAL_MAX/ddN.W )
      betac = DIR_TOTAL_MAX/ddN.W;
    if ( verbose>=1 && betain!=betac ) {
      yap_message("beta changed from %lf to %lf due to Dirichlet constrains\n",
		  betain, betac);
      ddP.betatot = betac*ddN.W;
    }
  }
  if ( ddP.window>0 ) {
    if ( ddP.window>=ddN.DT )
      ddP.window = 0;
    ddP.window_left = 0;
    ddP.window_right = ddP.window;
  }
  if ( ddP.Tinit==0 )
    ddP.Tinit = ddN.T;
  if ( ddP.NGalpha!=NULL) {
    ddT[ParNGAlpha].ptr = ddP.NGalpha;
  }
  if ( ddP.NGbeta!=NULL) {
    ddT[ParNGBeta].ptr = ddP.NGbeta;
  }
 }

/*
 *   initialising or ddP.betatot is changed
 *      file = where to grab beta
 *      resstem = only set if beta to be saved as ".beta"
 *
 *   when .betac in use, reset .betac and .betapr[] based in .betatot
 *   otherwise if .betapr[] in use, renormalise to get .betatot
 */
void pctl_fixbeta(char *file, char *resstem) {
  int c;
  if ( ddP.PYbeta!=H_HPDD && !ddP.betapr ) 
    ddP.betapr = dvec(ddN.W);
  if ( ddP.PYbeta==H_None && file==NULL && resstem )
    ddP.betac = 1;
  if ( ddP.PYbeta!=H_HPDD && file ) {
    /*
     *   only read on initialisation
     */
    ddP.betac = 0;
    if ( strcmp(file,"data")==0 ) {
      /*
       *  initialise according to occurrence in data
       *  but with Laplace smoothing
       */
      for (c=0; c<ddN.W; c++)
        ddP.betapr[c] = 0.5;
      for (c=0; c<ddN.NT; c++) 
        ddP.betapr[ddD.w[c]]++;
    } else if ( strcmp(file,"uniform")==0 ) {
      /*
       *  initialise uniform
       */
      for (c=0; c<ddN.W; c++)
        ddP.betapr[c] = 1.0;
    } else {
      /*
       *  initialise from file
       */
      read_dvec(file,ddN.W,ddP.betapr);
    }
    if ( resstem ) {
      /*
       *  however, record the result used
       */
      char *fname;
      fname = yap_makename(resstem,".beta");
      write_dvec(fname,ddN.W,ddP.betapr);
      free(fname);
    }
  }
  if ( ddP.betac!=0 && ddP.PYbeta==H_None ) {
    /*
     *  .betac  and .betapr set from .betatot
     */
    assert(ddP.betatot>0);
    ddP.betac = ddP.betatot/ddN.W;
    assert(ddP.betapr>0);
    for (c=0; c<ddN.W; c++)
      ddP.betapr[c] = ddP.betac;
  } else if ( ddP.betapr ) {
    /*
     *  renormalise .betapr[] so it adds to .betatot
     */
    double lastbeta = 0;
    for (c=0; c<ddN.W; c++)
      lastbeta += ddP.betapr[c];
    for (c=0; c<ddN.W; c++)
      ddP.betapr[c] *= ddP.betatot/lastbeta;
  }
}

/*
 *   initialising or ddP.alphatot is changed
 *      file = where to grab alpha
 *      resstem = only set if alpha to be saved as ".alpha"
 *
 *   when .alphac in use, reset .alphac and .alphapr[] based in .alphatot
 *   otherwise if .alphapr[] in use, renormalise to get .alphatot
 */
void pctl_fixalpha(char *file, char *resstem) {
  int c;
  assert(!ddP.NGalpha);
  if ( ddP.PYalpha!=H_HPDD && !ddP.alphapr ) 
    ddP.alphapr = dvec(ddN.T);
  if ( ddP.PYalpha==H_None && file==NULL && resstem )
    ddP.alphac = 1;
  if ( ddP.PYalpha!=H_HPDD && file ) {
    /*
     *   only read on initialisation
     */
    ddP.alphac = 0;
    if ( strcmp(file,"uniform")==0 ) {
      /*
       *  initialise uniform
       */
      for (c=0; c<ddN.T; c++)
        ddP.alphapr[c] = 1.0;
    } else {
      /*
       *  initialise from file
       */
      read_dvec(file,ddN.T,ddP.alphapr);
    }
    if ( resstem ) {
      /*
       *  however, record the result used
       */
      char *fname;
      fname = yap_makename(resstem,".alpha");
      write_dvec(fname,ddN.T,ddP.alphapr);
      free(fname);
    }
  }
  if ( ddP.alphac!=0 && ddP.PYalpha==H_None ) {
    /*
     *  .alphac  and .alphapr set from .alphatot
     */
    assert(ddP.alphatot>0);
    ddP.alphac = ddP.alphatot/ddN.T;
    assert(ddP.alphapr>0);
    for (c=0; c<ddN.T; c++)
      ddP.alphapr[c] = ddP.alphac;
  } else if ( ddP.alphapr ) {
    /*
     *  renormalise .alphapr[] so it adds to .alphatot
     */
    double lastalpha = 0;
    for (c=0; c<ddN.T; c++)
      lastalpha += ddP.alphapr[c];
    for (c=0; c<ddN.T; c++)
      ddP.alphapr[c] *= ddP.alphatot/lastalpha;
  }
}

/*
 *    called before the general dimensions have been set;
 *    as some parameters need to know dimensions to be set
 *    further fixes done by pctl_dims() and pctl_fixbeta()
 */
void pctl_fix(int ITER, int loadphi) {
  ddP.betapr = NULL;
  if ( ddP.ad==0 ) {
    ddT[ParAD].fix = 1;
  }
  if ( ddP.bdk==NULL) {
    ddT[ParBDK].fix = 1;   
    ddT[ParAD].fix = 1;
  }
  if ( ddP.PYalpha==H_NG ||  ddP.NGalpha ) {
    ddT[ParAlpha].fix = 1;
  } else {
    ddT[ParNGAlpha].fix = 1;
    ddT[ParNGBeta].fix = 1;
  }
  if ( ddP.PYalpha==H_None || ddP.PYalpha==H_NG || ddP.NGalpha ) {
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
    if ( ddP.a0==0 )
      ddT[ParA0].fix = 1;    
    if ( ddP.PYalpha==H_PDP ) {
      ddT[ParA0].fix = 1;    
      ddT[ParB0].fix = 1;   
    }
    if ( ddP.PYalpha==H_HDP||ddP.PYalpha==H_PDP ) 
      /*
       *    in this case .alphapr[] must be a prob vector
       */
      ddP.alphatot = 1;
    ddP.alphac = 0;
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
    }
    if ( ddP.PYbeta==H_HDP||ddP.PYbeta==H_PDP ) 
      /*
       *    in this case .betapr[] must be a prob vector
       */
      ddP.betatot = 1;
    ddP.betac = 0;
    ddT[ParBeta].fix = 1;
  }
  if ( loadphi ) {
    /*
     *   PYtheta and beta are not used!!
     */
    ddT[ParAW].fix = 1;
    ddT[ParAW0].fix = 1;
    ddT[ParBW].fix = 1;
    ddT[ParBW0].fix = 1;
    ddT[ParBeta].fix = 1;
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
  if ( ddP.mltiter>0 && ddP.hold_all==0 ) {
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

void pctl_report() {
  yap_message("PYbeta  = %d\n", (int)ddP.PYbeta);
  if ( ddP.betapr && ddP.betac==0 ) 
    yap_message("# beta proportions read from file\n");
  if ( ddP.PYbeta ) {
    yap_message("aw     = %lf\n", ddP.awpar);
    yap_message("bw     = %lf\n", ddP.bwpar);
    if ( ddP.PYbeta!=H_PDP ) {
      yap_message("aw0     = %lf\n", ddP.aw0);
      yap_message("bw0     = %lf\n", ddP.bw0);
    }
  } else {
    if ( ddP.betatot>0 )
      yap_message("betatot  = %lf # total over W=%d words\n", ddP.betatot, ddN.W);
  }
  yap_message("PYalpha  = %d\n", (int)ddP.PYalpha);
  if ( ddP.alphapr && ddP.alphac==0 ) 
    yap_message("# alpha proportions read from file\n");
  if ( ddP.PYalpha ) {
    yap_message("a     = %lf\n", ddP.apar);
    yap_message("b     = %lf\n", ddP.bpar);
    if ( ddP.PYalpha!=H_PDP ) {
      yap_message("a0     = %lf\n", ddP.a0);
      yap_message("b0     = %lf\n", ddP.b0);
    }
  } else if ( !ddP.NGalpha ) {
    yap_message("alphatot  = %lf # total over topics\n", 
		ddP.alphatot);
  }
  if ( ddP.bdk!=NULL ) {
    int t;
    yap_message("bdk =");
    for (t=0; t<ddN.T; t++) 
      yap_message(" %5lf", ddP.bdk[t]);
    yap_message("\n");
  }
  if ( ddP.NGalpha ) {
    int t;
    yap_message("NGalpha =");
    for (t=0; t<ddN.T; t++) 
      yap_message(" %5lf", ddP.NGalpha[t]);
    yap_message("\n");
    yap_message("NGbeta =");
    for (t=0; t<ddN.T; t++) 
      yap_message(" %5lf", ddP.NGbeta[t]);
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
 *   generate parameter corresponding to index;
 *   'iter' is cycle used to find which are active;
 *   and return in *par and *k 
 *        note bdk/NGalpha/NGbeta has ddP.kbatch values when used
 *        all other pars are single valued
 *   return 1 if found OK, else return 0 if no more
 */
int pctl_par_iter(int index, int iter, enum ParType *par, int *k) {
  enum ParType p;
  for (p=ParA; p<=ParBeta; p++) {
    if (  !ddT[p].fix && ddT[p].ptr
          && iter>ddT[p].start
	  && iter%ddT[p].cycles==ddT[p].offset ) {
      if ( p==ParBDK || p==ParNGAlpha || p==ParNGBeta ) {
        if ( index<ddP.kbatch) {
          *par = p;
          *k = (iter*ddP.kbatch/ddT[p].cycles+index)%ddN.T;
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
  int index;
  struct pst_data pd;
#ifdef H_THREADS
  int p;
  pthread_t thread[procs];
#endif
  
  /*
   *  first, create docstats if needed
   */
  ddP.docstats = NULL;
  //   why 100000?
  for (index=0; index<100000; index++) {
    int k;
    enum ParType par;
    int try = pctl_par_iter(index, iter, &par, &k);
    if ( try && par==ParBDK ) {
      ddP.docstats = dmi_bstore(&ddM);
    }
    if ( !try )
      break;
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
  if ( ddP.hold_all )
    return 1;
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
  yap_message("Sampling pars:");
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
  if ( ddP.NGbeta ) {
    /*
     *   give NGbeta an average of 1
     */
    int t;
    double tot = 0;
    for (t=0; t<ddN.T; t++)
      tot += ddP.NGbeta[t];
    tot /= ddN.T;
    for (t=0; t<ddN.T; t++)
      ddP.NGbeta[t] /= tot;
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
  fprintf(fp, "PYbeta  = %d\n", (int)ddP.PYbeta);
  if ( ddP.PYbeta ) {
    printpar(fp,ParAW); printpar(fp,ParBW);
    if ( ddP.PYbeta !=H_PDP ) {
      printpar(fp,ParAW0); printpar(fp,ParBW0);
    }
  } else {
    fprintf(fp, "#  %s is the total over W=%d words\n",ddT[ParBeta].name,ddN.W);
    printpar(fp,ParBeta);
  }
  fprintf(fp, "PYalpha  = %d\n", (int)ddP.PYalpha);
  if ( ddP.PYalpha ) {
    printpar(fp,ParA); printpar(fp,ParB);
    if ( ddP.PYalpha!=H_PDP ) {
      printpar(fp,ParA0); printpar(fp,ParB0);
    }
  } else if ( ddP.NGalpha ) {
    int t;
    if ( !ddT[ParNGAlpha].fix ) 
      fprintf(fp, "#  %s was sampled every %d major cycles in batches of %d\n", 
	      ddT[ParNGAlpha].name, ddT[ParNGAlpha].cycles, ddP.kbatch);
    fprintf(fp, "NGalpha =");
    for (t=0; t<ddN.T; t++) 
      fprintf(fp, " %5lf", ddT[ParNGAlpha].ptr[t]);
    fprintf(fp, "\n");
    if ( !ddT[ParNGBeta].fix ) 
      fprintf(fp, "#  %s was sampled every %d major cycles in batches of %d\n", 
	      ddT[ParNGBeta].name, ddT[ParNGBeta].cycles, ddP.kbatch);
    fprintf(fp, "NGbeta =");
    for (t=0; t<ddN.T; t++) 
      fprintf(fp, " %5lf", ddT[ParNGBeta].ptr[t]);
    fprintf(fp, "\n");
  } else {
    fprintf(fp, "#  %s is the total over topics\n",ddT[ParAlpha].name);
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
  if ( ddP.NGbeta!=NULL )
    free(ddP.NGbeta);
  if ( ddP.NGalpha!=NULL )
    free(ddP.NGalpha);
  if ( ddP.bdk!=NULL )
    free(ddP.bdk);
  if ( ddP.phi ) {
    free(ddP.phi[0]); free(ddP.phi);
  }
  if ( ddP.theta ) {
    free(ddP.theta[0]); free(ddP.theta);
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
  if ( ddP.alphapr )
    free(ddP.alphapr);
}
