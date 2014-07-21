/*
 * Utilities
 * Copyright (C) 2014 Wray Buntine
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@monash.edu     
 */

#include <stdio.h>  
#include <unistd.h>   
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "yap.h"
#include "fvec.h"

void fv_copy(float *v1, float *v2, int N) {
  memcpy(v1, v2, sizeof(v1[0])*N);
}
/*
 *  Hellinger distance fr prob vecs
 *    \sum_i (sqrt(p[i])-sqrt(q[i]))^2 
 *     = 2 - 2\sum_i sqrt(p[i]q[i])
 */
double fv_helldistunif(float *vp, int N) {
  double dist = 0;
  int i;
  for (i=0; i<N; i++ ) 
    dist += sqrt(vp[i]/N);
  return 2.0 * (1-dist);
}
double fv_helldist(float *vp, float *vp2, int N) {
  double dist = 0;
  int i;
  for (i=0; i<N; i++ ) 
    dist += sqrt(vp[i]*vp2[i]);
  return 2.0 * (1-dist);
}
/*
 *   noramlise too, just in case
 */
double fv_entropy(float *vp, int N) {
  double ent = 0;
  double tot = 0;
  int i;
  if ( !vp ) 
    return INFINITY;
  for (i=0; i<N; i++ ) {
    double p = vp[i];
    if ( N*p>0.00001 ) {
      tot += p;
      ent -= p * log(p);
    }
  }
  if ( tot<=0 )
    return INFINITY;
  return (ent+log(tot))/tot;
}
double fv_bound(float *vp, int N, float alpha) {
  int i;
  double count = 0;
  if ( !vp ) 
    return 0;
  for (i=0; i<N; i++ ) {
    if ( vp[i]>alpha )
      count ++;
  }
  return count;
}
double fv_avestrlen(float *vp, char **str, int N) {
  double asl = 0;
  int i;
  if ( !str || !vp )
    return 0.0;
  for (i=0; i<N; i++ ) {
    if ( str[i] )
      asl += strlen(str[i]) * vp[i];
  }
  return asl;
}
