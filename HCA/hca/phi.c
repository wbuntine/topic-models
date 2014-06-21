/*
 * Computing, updating and saving topicXword/phi estimates
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
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "yap.h"
#include "util.h"
#include "hca.h"
#include "data.h"
#include "stats.h"
#include "diag.h"
#include "check.h"


static char *phi_file = NULL;
static char *alpha_file = NULL;

void phi_init(char *resstem) {
  phi_file = yap_makename(resstem,".phi");
  ddG.phi_cnt = 0;
  if ( ! ddP.memory ) {
    ddS.phi = fmat(ddN.T,ddN.W);
  }
}
void alpha_init(char *resstem) {
  alpha_file = yap_makename(resstem,".alpha");
  ddG.alpha_cnt = 0;
  if ( ! ddP.memory ) {
    ddS.alpha = fvec(ddN.T);
  }
}

void phi_free() {
  ddG.phi_cnt = 0;
  if ( phi_file ) {
    free(phi_file);
    phi_file = NULL;
  }
  if ( ddS.phi ) {
    free(ddS.phi[0]); free(ddS.phi);
    ddS.phi = NULL;
  }
}
void alpha_free() {
  ddG.alpha_cnt = 0;
  if ( alpha_file ) {
    free(alpha_file);
    alpha_file = NULL;
  }
  if ( ddS.alpha ) {
    free(ddS.alpha);
    ddS.alpha = NULL;
  }
}

/*
 *    only need to save if *not* storing in file
 */
void phi_save() {
  if ( !ddP.memory ) {
    uint32_t cnt, dim;
    FILE *fpout=NULL;
    fpout = fopen(phi_file,"wb");
    if ( !fpout ) 
      yap_sysquit("Cannot open file '%s' for write in phi_save()\n", 
		  phi_file);
    dim = ddN.W;
    if ( fwrite(&dim, sizeof(dim), 1, fpout) !=1 )
      yap_sysquit("Cannot write dim to '%s' in phi_save()\n", phi_file);
    dim = ddN.T;
    if ( fwrite(&dim, sizeof(dim), 1, fpout) !=1 )
      yap_sysquit("Cannot write dim to '%s' in phi_save()\n", phi_file);
    cnt = ddG.phi_cnt+1;
    if ( fwrite(&cnt, sizeof(cnt), 1, fpout) !=1 )
      yap_sysquit("Cannot write count to '%s' in phi_save()\n", phi_file);
    if ( fwrite(ddS.phi[0], sizeof(ddS.phi[0][0]), ddN.W*ddN.T, fpout) 
	 !=ddN.W*ddN.T )
      yap_sysquit("Cannot write matrix to '%s' in phi_save()\n", phi_file);
    fclose(fpout);
#ifndef NDEBUG
    {
      int w,t;
      for (t=0; t<ddN.T; t++)
	for (w=0; w<ddN.W; w++) 
	  if ( ddS.phi[t][w]<=0 )
	    yap_message(" saving ddS.phi[%d][%d]<=0\n", t, w);
    } 
#endif
  }
}

void alpha_save() {
  write_fvec(alpha_file,ddN.T,ddS.alpha);
  if ( verbose>1 ) {
    int t;
    yap_message("Saving ddS.alpha:");
    for (t=0; t<ddN.T; t++)
	yap_message(" %f", ddS.alpha[t]);
    yap_message("\n");
  }
}

double phi_entropy(int k) {
  double ent = 0;
  int i;
  assert(ddP.phi || ddS.phi);
  if ( ddP.phi ) {
    for (i=0; i<ddN.W; i++ ) {
      double p = ddP.phi[k][i];
      if ( p>0.00001 ) 
	ent -= p * log(p);
    }
  } else {
   for (i=0; i<ddN.W; i++ ) {
      double p = ddS.phi[k][i];
      if ( p>0.00001 ) 
	ent -= p * log(p);
   }
  }
  return ent;
}
double alpha_entropy() {
  double ent = 0;
  int k;
  assert(ddP.fixalpha || ddS.alpha);
  if ( ddP.fixalpha ) {
    for (k=0; k<ddN.T; k++) {
      double p = ddP.fixalpha[k];
      if ( p>0.00001 ) 
	ent -= p * log(p);
    }
  } else {
    for (k=0; k<ddN.T; k++) {
      double p = ddS.alpha[k];
      if ( p>0.00001 ) 
	ent -= p * log(p);
    }
  }
  return ent;
}
 
void phi_load(char *resstem) {
  uint32_t cnt, dim;
  int totgot, got;
  FILE *fp=NULL;
  phi_file = yap_makename(resstem,".phi");
  ddP.phi = fmat(ddN.T,ddN.W);
  fp = fopen(phi_file,"rb");
  if ( !fp ) 
    yap_sysquit("Cannot open file '%s' for read in phi_load()\n", 
		phi_file);
  if ( fread(&dim, sizeof(dim), 1, fp) !=1 )
    yap_sysquit("Cannot read dim to '%s' in phi_load()\n", phi_file);
  if ( dim != ddN.W )
    yap_quit("Bad word dim in phi_load()\n");
  if ( fread(&dim, sizeof(dim), 1, fp) !=1 )
    yap_sysquit("Cannot read dim to '%s' in phi_load()\n", phi_file);
  if ( dim != ddN.T ) 
    yap_quit("Bad topic dim in phi_load()\n");
  if ( fread(&cnt, sizeof(cnt), 1, fp) !=1 )
    yap_sysquit("Cannot read count to '%s' in phi_load()\n", phi_file);
  ddG.phi_cnt = cnt;
  got = 0;
  totgot = 0;
  while ( (got=(int)fread(&ddP.phi[0][got], sizeof(ddS.phi[0][0]), ddN.W*ddN.T-totgot, fp)) 
       !=ddN.W*ddN.T-totgot ) {
    yap_message("load_phi(%d,%d): got %d in %d\n", (int)ddN.W, (int)ddN.T, got, totgot);
    if ( got==0 || ferror(fp) || feof(fp) )
       yap_sysquit("Cannot read matrix from '%s' in phi_load()\n", phi_file);
    totgot += got;
  }
  totgot += got;
  assert(totgot==ddN.W*ddN.T);
  fclose(fp);
  if ( verbose )
     yap_message("Read matrix from '%s' in phi_load()\n", phi_file);
  free(phi_file);
  phi_file = NULL;
#ifndef NDEBUG
  if ( verbose>0 ) {
    /*  check how well normalises */
    double tot;
    int w, t;
    int bad = 0;
    for (t=0; t<ddN.T; t++) {
      tot = 0;
      for (w=0; w<ddN.W; w++)
	tot += ddP.phi[t][w];
      if ( fabs(tot-1.0)>1e-5 && bad<5 ) {
	yap_message("phi_load():  ddP.phi[%d][*] sums to %lf\n", t, tot);
	bad++;
	break;
      }
    }
  }
#endif
}

void phi_update() {
  int t, w;
  static FILE *fpin=NULL;
  static FILE *fpout=NULL;
  float *vec = NULL;

  if ( ddP.memory ) {
    if ( ddG.phi_cnt>0 ) {
      char *fname=yap_makename(phi_file,".tmp");
      if ( rename(phi_file,fname) )
	yap_sysquit("Cannot rename '%s' in phi_update()\n",phi_file);
      fpin = fopen(fname,"rb");
      if ( !fpin )
	yap_sysquit("Cannot open '%s' for read in phi_update()\n",fname);
      free(fname);
    } else {
      /*  delete any existing file */
      unlink(phi_file);
    }
    vec = fvec(ddN.W);
    fpout = fopen(phi_file,"wb");
    if ( !fpout ) 
      yap_sysquit("Cannot open file '%s' for write in phi_update()\n",
		  phi_file);
  } 
  /*
   *   files now open, but fpin==NULL if none there
   */
  if ( ddP.memory ) {
    /*
     *  read/write headers
     */
    uint32_t dim, cnt;
    if ( fpin ) {
      uint32_t dims[3];
      if ( fread(dims, sizeof(dims[0]), 3, fpin) !=3 )
	yap_sysquit("Cannot read from '%s.tmp' in phi_update()\n", phi_file);
      if ( dims[0]!=ddN.W ||  dims[1]!=ddN.T )
	yap_quit("Bad dimensions in phi_update()\n");
    }
    dim = ddN.W;
    if ( fwrite(&dim, sizeof(dim), 1, fpout) !=1 )
      yap_sysquit("Cannot write dim to '%s' in phi_update()\n", phi_file);
    dim = ddN.T;
    if ( fwrite(&dim, sizeof(dim), 1, fpout) !=1 )
      yap_sysquit("Cannot write dim to '%s' in phi_update()\n", phi_file);
    cnt = ddG.phi_cnt+1;
    if ( fwrite(&cnt, sizeof(cnt), 1, fpout) !=1 )
      yap_sysquit("Cannot write count to '%s' in phi_update()\n", phi_file);
  }

  /*
   *    now compute
   */
  for (t=0; t<ddN.T; t++) {
    double totvec = 0;
    if ( ddP.memory ) {
      if ( fpin ) {
	if ( fread(vec, sizeof(vec[0]), ddN.W, fpin) !=ddN.W )
	  yap_sysquit("Cannot read vector from '%s.tmp' in phi_update()\n", 
		      phi_file);
      }
    } else {
      vec = ddS.phi[t];
    }
    if ( ddG.phi_cnt==0 ) {
      for (w=0; w<ddN.W; w++) vec[w] = 0;
    }
    
    for (w=0; w<ddN.W; w++) {
      double val = wordprob(w, t);
      totvec += val;
      if ( val<=0 )
	yap_message("phi_update for phi[%d][%d] zero\n",t,w);
      vec[w] = (ddG.phi_cnt*vec[w] + val) / (ddG.phi_cnt+1);
    }
#ifndef NDEBUG
    if ( fabs(totvec-1.0)>0.02 )  {
      double bwptot = 0;
      assert(ddP.PYbeta==H_HPDD) ;
      yap_message("Word vector for topic %d doesn't normalise (%g)\n",
		  t, totvec);
      for (w=0; w<ddN.W; w++) 
	bwptot += betabasewordprob(w);
      yap_message("  \\sum_w betabasewordprob(w) = %f\n", bwptot);
      check_Tw();
    }
#endif
    if ( ddP.memory ) {
      if ( fwrite(vec, sizeof(vec[0]), ddN.W, fpout) !=ddN.W )
	yap_sysquit("Cannot write vector to '%s' in phi_update()\n", 
		    phi_file);
    } else {
      vec = NULL;
    }
  }

  /*
   *  close up
   */
  ddG.phi_cnt++;
  if ( fpin ) {
    char *fname=yap_makename(phi_file,".tmp");
    fclose(fpin);
    if ( unlink(fname) )
      yap_sysquit("Cannot unlink '%s' in phi_update()\n",fname);
    free(fname);
  }
  if ( fpout )
    fclose(fpout);
  if ( vec )
    free(vec);
}

void alpha_update() {
  int t;
  double totvec;
  
  /*
   *    now compute
   */
  totvec = 0;
  for (t=0; t<ddN.T; t++) {
    double val = alphabasetopicprob(t);
    totvec += val;
    if ( val<=0 )
      yap_message("alpha_update for alpha[%d] zero\n",t);
    ddS.alpha[t] = (ddG.alpha_cnt*ddS.alpha[t] + val) / (ddG.alpha_cnt+1);
  }
#ifndef NDEBUG
  if ( fabs(totvec-1.0)>0.02 )  {
    yap_message("Base probability for topics doesn't normalise (%g)\n",
		totvec);
  }
#endif
  ddG.alpha_cnt++;
}
