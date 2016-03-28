/*
 * Auxiliary topic routines for predictiing top k
 * Copyright (C) 2016 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include "hca.h"
#include "yap.h" 

static float *wvec;
static int pcompar(const void *a, const void *b) {
  float na = wvec[*(uint32_t*)a];
  float nb = wvec[*(uint32_t*)b];
  if ( na<nb ) return 1;
  if ( na>nb ) return -1;
  return 0;
}

void predict_topk(char *resstem, int topk) {
  char *fname;
  FILE *fp;
  
  uint32_t *uvec;  /* for sorting */
  int i,k,w,l;
  fname = yap_makename(resstem,".topk");
  fp = fopen(fname,"w");
  if ( !fp ) {
    yap_sysquit("Cannot open '%s' for write\n", fname);
  }
  wvec = fvec(ddN.W);
  uvec = u32vec(ddN.W);
  for (i=0; i<ddN.DT; i++) {
    //  build word probs
    for (w=0; w<ddN.W; w++) 
      wvec[w] = 0;
    for (k=0; k<ddN.T; k++) {
      float ww = ddP.theta[i][k];
      for (w=0; w<ddN.W; w++)
	wvec[w] += ww*ddP.phi[k][w];
    }
    //  remove words already seen
    for (l=ddD.NdTcum[i]; l<ddD.NdTcum[i+1]; l++) 
      wvec[ddD.w[l]] = 0;
    //  sort
    for (w=0; w<ddN.W; w++) 
      uvec[w] = w;
    qsort(uvec, ddN.W, sizeof(*wvec), pcompar);
    //   print
    for (w=0; w<topk; w++)
      fprintf(fp, " %u", uvec[w]);
      // fprintf(fp, " %u(%g)", uvec[w], wvec[uvec[w]]);
    fprintf(fp, "\n");
  }
  free(wvec);
  free(uvec);
  fclose(fp);
  free(fname);
}
