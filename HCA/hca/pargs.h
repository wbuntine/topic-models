/*
 * arguments for parallel call to Gibbs samplers
 * Copyright (C) 2014 Swapnil Mishra
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Swapnil Mishra (swapnil.mishra@anu.edu.au)
 *
 *
 */
#ifndef __PARGS_H
#define __PARGS_H

#include "gibbs.h"


/*
 *   arguments for parallel call to Gibbs samplers
 */
typedef struct D_pargs_s {
  enum GibbsType fix;
  int Tmax;
  int window;
  double thislp;
  int thisNd;
  int dots;
  int processid;
  int procs;  
  double tot_time;
  int *doc;
} D_pargs_p;
#endif
