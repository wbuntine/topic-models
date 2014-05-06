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
#include "pctl.h"
#include "stats.h"
#include "srng.h"

#define MAXM 1000
/* 
 *   when defined, stops introducing new topics into a 
 *   document after the first ... ???
 */
// #define GIBBS_ZEROD

/* 
 *   Allow experimental stuff ... only for developers since
 *   they mostly don't work
 */
// #define EXPERIMENTAL

/*
 *   Switch on to allow threading
 *   if off some vestiges remain but wont call threads
 */
#define H_THREADS

/*
 *   define this to set up the cache in alphabasetopicprob()
 *   doesn't work with threads
 */
#ifndef H_THREADS
// #define CACHE_ABTP
#endif

/*
 *   when defined does tracking of changes to a single Nwt
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

/*
 *    type of prior node for PDP
 *      - none
 *      - hierarchical DP
 *      - hierarchical PDD (GEM-like)
 *      - constant
 */
enum PDPType { H_None=0, H_HDP, H_HPDD, H_PDP };

double likelihood();

double lp_test_Pred(char *resstem);
#ifdef EXPERIMENTAL
double lp_test_LRS();
#endif
double lp_test_ML(int procs, enum GibbsType fix);

void query_read(char *fname);
void gibbs_query(int K, char *qname, int dots, int this_qpart, int qparts);

void print_maxz(char *fname);

//==================================================
// global variables
//==================================================

extern rngp_t rngp;
extern int verbose;

#endif
