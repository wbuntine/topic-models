/*
 * Sparse word recording 
 * Copyright (C) 2013 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *
 *     Two cases:
 *           some words are recorded, so we have special bitmaps
 *                to support checking
 *           all words are being recorded, in which case we need a 
 *               short-circuit of the above
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include "hca.h"
#include "stats.h"
#include "util.h" 
#include "yap.h" 
#include "diag.h" 

/*
 *    report sparsity for words in map
 *    if fp is NULL do it for all words!
 */
void sparsemap_init(FILE *fp, int procs) {
  int w, n=0;
  /*
   *    set up words
   */
  if ( fp ) {
    int dim = 10;
    ddG.words = malloc(sizeof(*ddG.words)*dim);
    while ( fscanf(fp," %d", &w)==1 ) {
      if ( n>=dim ) {
	dim+=10;
	ddG.words = realloc(ddG.words, dim*sizeof(*ddG.words));
      if ( !ddG.words )
	yap_quit("Out of memory in sparsemap_init()\n");
      }
      ddG.words[n++] = w;
    }
    //  ddG.words = realloc(ddG.words, n);
    ddG.n_words = n;
  } else {
    ddG.words = NULL;
    ddG.n_words = ddN.W;
  }
 /*
   *   set up data structures
   */
  ddG.code = malloc(sizeof(*ddG.code)*procs);
  if ( !ddG.code )
    yap_quit("Out of memory in sparsemap_init()\n");
  {
    int p;
    for (p=0; p<procs; p++)
      ddG.code[p] = fmat(ddG.n_words,ddN.T);
  }
  if ( fp ) {
    ddG.iscodeword = u32vec((ddN.W+31)/32);
    for (n=0; n<ddG.n_words; n++) {
      ddG.iscodeword[ddG.words[n]/32U] |= (1U << (ddG.words[n]%32U));
    }
  } else {
    ddG.iscodeword = NULL;
  }
  ddG.docode = 0;
  ddG.didcode = 0;
}

void sparsemap_null() {
  ddG.n_words = 0;
  ddG.words = NULL;
  ddG.code = NULL;
  ddG.iscodeword = NULL;
  ddG.docode = 0;
}

void sparsemap_free() {
  if ( ddG.n_words>0 ) {
    ddG.n_words = 0;
    if ( ddG.words )
      free(ddG.words);
    free(ddG.code[0][0]);free(ddG.code[0]);free(ddG.code);
    if ( ddG.iscodeword )
      free(ddG.iscodeword);
    sparsemap_null();
  }
}

void sparsemap_report(char *resstem, double epsilon, int procs) {
  FILE *fp;
  int n, k;
  char *fname = yap_makename(resstem,".smap");
  fp = fopen(fname,"w");
  if ( !fp )
    yap_sysquit("Cannot open sparse map output file '%s'\n", fname);
  for (n=0; n<ddG.n_words; n++) {
    double ent = 0;
    double tot = 0;
    uint32_t w = n;
    if ( ddG.words ) 
	w = ddG.words[n];
    fprintf(fp,"%s(%d):", ddN.tokens?ddN.tokens[w]:"--", w);
    if ( procs>1 ) {
      /*  accummulate into ddG.code[0]  */
      int p;
      for (k=0; k<ddN.T; k++) 
	for (p=1; p<procs; p++) 
	  ddG.code[0][n][k] += ddG.code[p][n][k];
    }
    for (k=0; k<ddN.T; k++) {
      double val = ddG.code[0][n][k]/ddG.didcode;
      if ( val>epsilon ) 
	fprintf(fp," %d/%0.1f", k, val);
      tot += val;
    }
    for (k=0; k<ddN.T; k++) {
      double val = ddG.code[0][n][k]/ddG.didcode/tot;
      if ( val>0.00001 ) 
	ent -= val * log(val);
    }
    fprintf(fp," perp=%lf\n", exp(ent));
  }
  fclose(fp);
  free(fname);
}

int sparsemap_word(uint32_t w) {
  int n;
  /*  short circuit if all words recorded */
  if ( !ddG.iscodeword )
    return w;
  /*  otherwise, lookup */
  assert(G_isword(w));
  for (n=0; n<ddG.n_words; n++)
    if ( ddG.words[n]==w )
      break;
  return n;
}
