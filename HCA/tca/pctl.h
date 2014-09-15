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
#define B0PAR 10
#define AWPAR 0.5
#define BWPAR 5000
#define ACYCLES 11      //  by default update a's
#define BCYCLES 3      //  by default update b's

/*
 *    set to record stats of indicator values each cycle
 *    to RESSTEM.istats
 */
// #define IND_STATS

/*
 *  hyperparameters and test parameters
 */
typedef struct D_pars_s {
  /*
   *  parameters, only used for testing
   */
  float ***phi;
  float **mu;
  /*
   *    hyperparameters
   */
  double a_mu;             //  
  double *b_mu;            //   indexed by e
  double b_mu0;            //  
  double a_phi0;            //  for phi at epoch 0
  double a_phi1;            //  for phi at others
  double **b_phi;          //   indexed by e, t
  double b_phi0;           //   one for the prior \phi_0
  double a_theta;          // for mu->Document step (mid)
  double b_theta;          // for mu->Document step (mid)
  double a_burst, b_burst;     //  burstiness

  /*
   *  test and report controls
   */
  int phiiter, phiburn;     // burnin and iterations for topicXword probs
  int muiter, muburn;       // burnin and iterations for topic prior probs
  int mltiter, mltburn;     //  burnin and iterations for ML testing
  char *cofile;             //  set if want to do PMI-based coherency test
  int progiter, progburn;   //  progress reports
  int training;             //  suggested training set size
  char *teststem;           //  stem for the test data, only if different
  /*
   *  special control for sampling batches of pars
   */
  int kbatch;
  /*
   *  bound controlling how far back a table can go
   */
  int back;
  /*
   *   bounds for Stirling number tables
   */
  int maxN;
  int maxM;
  /*
   *  hold out method
   *     hold_every>0 then hold out every n-th in doc
   *     hold_dict>0 then hold out every n-th in dictionary
   *     hold_fraction>0 then hold out final part
   */
  double hold_fraction;
  int    hold_dict;
  int    hold_every;
#ifdef IND_STATS
  /*
   *   indicator stats
   */
  uint32_t ***doc_ind_stats;
  uint32_t ***word_ind_stats;
#endif
} D_pars_t;

#define PCTL_BURSTY()          (ddP.b_burst>0)
#define PCTL_NOMU()  (ddP.b_mu[0]>PYP_CONC_MAX*100)

/*
 *  hyperparameters control
 */
enum ParType { ParNone=0,
	       ParAM, ParBM0, ParBM1, /*   for a_mu, b_mu[] */
	       ParB0M,                /*   for b_mu0 */
	       ParAP0, ParAP1, ParBP0, ParBP1, /*   for a_phi, b_phi[][] */
	       ParB0P,                /*   for b_phi0 */
	       ParAT, ParBT,          /*   for a_theta, b_theta */
	       ParAB, ParBB           /*   for a_burst, b_burst */
             };

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

extern D_pars_t ddP;
extern D_pctl_t ddT[];
enum ParType findpar(char *name);
void pctl_init();
void pctl_read(char *resstem, char *buf);
void pctl_fix(int ITER);
void pctl_report();
void pctl_sample(int iter, int procs);
void pctl_update(int iter);
void pctl_print(FILE *fp);
void pctl_samplereport();
int pctl_hold(int i);
int pctl_training(int D);
void pctl_free();

double pctl_gammaprior(double x);
void pctl_loadphi(char *resstem);
void pctl_loadmu(char *resstem);

#endif
