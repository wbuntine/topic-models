/*
 * Probability utilities
 * Copyright (C) 2010-2011 Wray Buntine
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
#include <time.h>

#include "yap.h"
#include "util.h"
#include "data.h"
#include "tca.h"
#include "stats.h"


/*
 *   no fancy statistics here, just the empirical
 *   topic counts smoothed with alpha=0.1;
 *   in fact its meaningless statistically, just a
 *   useful diagnostic
 */
void get_probs(double *vp) {
  int t, d, NWT=0;
  for (d=0; d<ddN.DT; d++)
    NWT += ddS.N_dT[d];
  for (t=0; t<ddN.T; t++) {
    int tot=0;
    for (d=0; d<ddN.DT; d++)
      tot += ddS.n_dt[d][t];
    vp[t] = (0.1+tot)/(0.1*ddN.T+NWT);
  }
}

void yap_probs() {
  int t;
  int empty = 0;
  double ent = 0;
  double *vp = malloc(sizeof(*vp)*ddN.T);
  get_probs(vp);
  yap_message("probs = ");
  for (t=0; t<ddN.T; t++) 
    if ( vp[t]>0 ) {
      yap_message(" %lf", vp[t]);
      ent -= vp[t]*log(vp[t]);
    } else {
      empty++;
      yap_message(" -");
    }
  yap_message("\nempty = %d, ent = %lf\n", empty, exp(ent));
  free(vp);
}

void print_probs(FILE *fp) {
  int t, num = 0;
  double *vp = malloc(sizeof(*vp)*ddN.T);
  get_probs(vp);
  fprintf(fp, "probs = ");
  for (t=0; t<ddN.T; t++) 
    if ( vp[t]>0 ) {
      fprintf(fp, " %lf", vp[t]);
      num++;
    } else
      fprintf(fp, " -");
  fprintf(fp, "\n");
  fprintf(fp, "# topics: %d\n", num);
  free(vp);
}

