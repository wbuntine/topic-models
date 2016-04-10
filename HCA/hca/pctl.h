/*
 * Hyperparameters, test parameters, and their control
 * Copyright (C) 2012-2014 Wray Buntine
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
#ifndef __PCTL_H
#define __PCTL_H

#include <stdint.h>
#include "gibbs.h"
#include "dirdim.h"

#define BETA 100
#define APAR 0.0
#define BPAR 10
#define A0PAR 0.0
#define B0PAR 10
#define AWPAR 0.5
#define BWPAR 100
#define AW0PAR 0.5
#define BW0PAR 10
#define ACYCLES 11      //  by default update a's
#define BCYCLES 3      //  by default update b's
#define DIRCYCLES 4     //  by default update alpha/beta
#define STARTCYCLES 1      //  when to start sampling hypers
/*
 *   pars for prior Gamma() for NG alpha
 */
#define NGASH 1.1
#define NGASC 1.0

/*
 *   pars for prior Beta()
 */
#define NGS0 1.0
#define NGS1 1.0

/*
 *  hyperparameters and test parameters
 */
typedef struct D_pars_s {
  /*
   *  sometimes we fix \phi, \theta or \alpha during testing (or training)
   */
  float **theta;
  float **phi;
  /*
   *    hyperparameters for alpha Dirichlet
   */
  double *alphapr;        // vector normalises to 1 when PYalpha!=H_None
  double alphac;          // individual constant, set from alphatot
  double alphatot;        // total of above or alphac*T

 /*
   *    hyperparameters for beta Dirichlet
   */
  double *betapr;         // vector normalises to 1 when PYbeta!=H_None
  double betac;           // individual constant, set from betatot
  double betatot;         // total of above or betac*W

  /*
   *  Pitman-Yor parameters
   */
  enum PDPType   PYalpha; // non-zero if using table indicators, 
  double apar, bpar;
  double a0, b0;          // PDD/PDP params for root
  enum PDPType   PYbeta;  // non-zero if using table indicators, 
  double awpar, bwpar;
  double aw0, bw0;        // PDD/PDP params for W root

  /*
   *  normalised gamma stuff, also uses flag in PYalpha
   */
  double *NGbeta;         // beta vector for H_NG
  double *NGalpha;        //  alpha vector
  double ngash, ngasc;    // shape and scale for Gamma prior of NG alpha
  double ngs0, ngs1;       // Beta priors for sparsity

  /*
   *  burstiness stuff
   */
  double ad;              // PDP params for doc
  double *bdk;            // version with seperate bd for each topic
  
  /*
   *   min. size of allowed doc
   */
  int mindocsize;
  /*
   *  test and report controls
   */
  int prditer, prdburn;     //  burnin and iterations for prediction tests
  int lrsiter, lrsburn;     //  burnin and iterations for LRS testing
  int mltiter, mltburn;     //  burnin and iterations for ML testing
  char *cofile;             //  set if want to do PMI-based coherency test
  int empirical;            //  =0 means '-l' records pars, else data
  int spiter, spburn;       //  burnin and iterations for sparsity testing
  int tprobiter, tprobburn;  //  burnin and iters for test docXtopic probs
  int probiter, probburn;   //  burnin and iters for train docXtopic probs
  int phiiter, phiburn;     // burnin and iterations for topicXword probs
  int alphaiter, alphaburn;     // burnin and iterations for topic prior probs
  int progiter, progburn;   //  progress reports
  int queryiter;            //  iterations for query
  int queryiter0;           //  iterations for query, first pass
  int memory;               //  higher value means conserve more memory
  int training;             //  suggested training set size
  char *teststem;           //  stem for the test data, only if different
  int mergeiter, mergeinit; //  when to start, number of iterations
  float mergemin;             //  ignore topics with this or less proportion
  float topcor;              // report topic correlation if GT
  int mergebest;            //  include best non-clashing merges at each round
  /*
   *     window control ... only work on this much data at once
   */
  int window;               //  size
  int window_incr;          //  change by this much each cycle
  int window_cycle;         //  cycle to begin moving
  int window_left;          //  bounds, is treated modulo
  int window_right;
  /*
   *  special control for sampling P.bdk[]
   */
  int kbatch;
  uint16_t **docstats;
  /*
   *  querying, multiple queries stored in single vector
   */
  uint32_t *qword;      /*  dictionary index for word */
  int16_t *query;      /*  map from dictionary index to first qword[] index */
  int16_t *qposn;      /*  index to qword{} where this query starts */
  uint32_t *qid;       /*  query number for word, base is 0 */

  int     n_excludetopic;
  int     *excludetopic;  /*  list of topics to exclude from prob */
  uint32_t *bits_et;      /*  boolean vector version */
  
  int n_query;         /*  count of queries */
  int n_words;         /*  count of words in all queries */
  /*
   *  incrementing topics ... a maximum topic count maintained
   */
  int Tinit;                //   starting number of topics allowed (0=max)
  int Tcycle;               //   cycles when allowed to change
  int Tinc;                 //   allowed increment
  int Tfree;                //   after this many cycles, drop constraints
  /*
   *  hold out method
   *     hold_every>0 then hold out every n-th in doc
   *     hold_dict>0 then hold out every n-th in dictionary
   *     hold_fraction>0 then hold out final part
   */
  double hold_fraction;
  int    hold_dict;
  int    hold_every;
  int    hold_all;
} D_pars_t;

#define PCTL_BURSTY()          (ddP.bdk!=NULL)

//  #define BWPAR0
#ifdef BWPAR0
#define ddP_bwpar(t)  (t==0?10000:ddP.bwpar)
#else
#define ddP_bwpar(t)  ddP.bwpar
#endif

/**************
 * 
 *  hyperparameters control:   all set up so hyperparameters can be
 *  sampled semi-automatically ... using C++ might be easier ;-)
 *
 */
enum ParType { ParNone=0, ParA, ParB, ParA0, ParB0, 
	       ParAW, ParBW, ParAW0, ParBW0, 
               ParAD, ParBDK, ParNGBeta, ParNGAlpha,
	       ParNGS0, ParNGS1, ParNGASH, ParNGASC,
	       ParAlpha, ParBeta };
typedef struct D_pctl_s {
  char *name;
  double *ptr;
  char fix;
  int start;
  int offset;
  int cycles;
  void (*sampler)(double *x);
  void (*samplerk)(double *x, int k);
} D_pctl_t;

#define Q_excludetopic(k) (ddP.bits_et[(k)/32U] & (1U << (((unsigned)k)%32U)))

extern D_pars_t ddP;
extern D_pctl_t ddT[];
enum ParType findpar(char *name);
/*
 *   initialisation sequence ... a mess
 */
/*  create blank pctl */
void pctl_init();  
void pctl_read(char *resstem);
/*  adjust pctl after restart/commandline  */
void pctl_fix(int ITER, int loadphi);
/*   have data+topic dims so adjust vectors */
void pctl_dims();
void pctl_fixalpha(char *file, char *resstem);
void pctl_fixbeta(char *file, char *resstem);

/*
 *   processing
 */
void pctl_report();
void pctl_sample(int iter, int procs);
void pctl_update(int iter);
void pctl_print(FILE *fp);
void pctl_samplereport();
int pctl_Tmax(int Tmax, int iter);
int pctl_hold(int i);
int pctl_training(int D);
void pctl_free();
void pctl_query(char *qname);

double pctl_gammaprior(double x);

double pctl_ng_alphapriorZ();
double pctl_ng_alphaprior(double x);
void pctl_ng_normbeta();

#endif
