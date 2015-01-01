/*
 * Module for reading and initialising term stats
 * Copyright (C) 2014 Wray Buntine
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@monash.edu)
 *
 */

#ifndef __TERMSTATS_H
#define __TERMSTATS_H

#include <stdint.h>

/*
 *    used to build global data about word unique usage for a doc
 */
typedef struct T_stats_s {
  /*
   *  basic dims and data stored/saved from elsewhere
   */
  int T;
  int DT;
  /*
   *  counts for terms
   *  terms are Kmin, Kmin+1, ..., Kmin+K-1
   */
  int K;  
  int Kmin;
  uint32_t **Nkt;
  char **tokens;
} T_stats_t;


void tstats_free(T_stats_t *ptr);
/*
 *      returns structure or NULL if none
 *
 *      all arguments come from the standard data structures
 */
T_stats_t *tstats_init(uint16_t *z, uint32_t *NdTcum, //  cumsum(NdT)
		       int T, int DT,  // dims
		       char *stem);

#endif
