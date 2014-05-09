/*
 * Basic definitions/types
 * Copyright (C) 2009-2014 Wray Buntine
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
#ifndef __HCA_H
#define __HCA_H

#include <unistd.h>
#include "lgamma.h" 
#include "util.h" 
#include "stable.h"
#include "srng.h" 
#include "pctl.h"
#include "stats.h"

#define MAXM 1000

/*
 *   Switch on to allow threading
 *   if off some vestiges remain but wont call threads
 *   NB. some asserts must be off with threading due to optimisation
 */
// #define H_THREADS

/*
 *   when defined does traing of changes to a single m_evt
 *   during sampling
 */
// #define TRACE_WT
#ifdef TRACE_WT
#define TR_W 4744
#define TR_T 7
#endif

/*
 *    used when printing words
 */
enum ScoreType { ST_count, ST_idf, ST_cost, ST_Q, ST_phi };

void cache_init() ;
void cache_free();

double likelihood();

double lp_test_ML(enum GibbsType fix);
void print_maxz(char *fname);

//==================================================
// global variables
//==================================================

extern rngp_t rngp;
extern int verbose;

/*
 *  Cache
 */
typedef struct D_cache_s {
  stable_t *a_mu;
  stable_t *a_phi;
  stable_t *a_theta;
  stable_t *a_burst;
} D_cache_t;

extern D_cache_t ddC;

#endif
