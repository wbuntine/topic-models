/*
 *   arguments for parallel function
 */
#ifndef __PARGS_H
#define __PARGS_H

#include "gibbs.h"


typedef struct D_pargs_s {
  enum GibbsType fix;
  int Tmax;
  double thislp;
  int thisNd;
  int dots;
  int processid;
  int procs;  
  double tot_time;
} D_pargs_p;
#endif
