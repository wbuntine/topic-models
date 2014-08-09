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
 *  If ddP.memory is set, then statistics kept
 *  in file in binary format.
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
static char *beta_file = NULL;

#define ESTIMATE_BETA (ddP.PYbeta && ddP.PYbeta!=H_PDP)

void phi_init(char *resstem) {
  phi_file = yap_makename(resstem,".phi");
  ddG.phi_cnt = 0;
  if ( ! ddP.memory ) {
    ddS.phi = fmat(ddN.T + (ESTIMATE_BETA?1:0),ddN.W);
  }
  if ( ESTIMATE_BETA ) 
    beta_file = yap_makename(resstem,".beta");
}
void alpha_init(char *resstem) {
  alpha_file = yap_makename(resstem,".alpha");
  ddG.alpha_cnt = 0;
  ddS.alpha = fvec(ddN.T);
}

void phi_free() {
  ddG.phi_cnt = 0;
  if ( phi_file ) {
    free(phi_file);
    phi_file = NULL;
  }
  if ( beta_file ) {
    free(beta_file);
    beta_file = NULL;
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
  free(ddS.alpha);
  ddS.alpha = NULL;
}

/*
 *    only need to save if *not* storing in file
 */
void phi_save() {
  if ( !ddP.memory ) {
    uint32_t cnt, dim, size;
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
    size = ddN.W*(ddN.T+(ESTIMATE_BETA?1:0));
    if ( fwrite(ddS.phi[0], sizeof(ddS.phi[0][0]), size, fpout) 
	 !=size )
      yap_sysquit("Cannot write matrix to '%s' in phi_save()\n", phi_file);
    fclose(fpout);
    if ( ESTIMATE_BETA ) {
      assert(beta_file);
      write_fvec(beta_file,ddN.W,ddS.phi[ddN.T]);
    }
#ifndef NDEBUG
    {
      int w,t;
      for (t=0; t<ddN.T; t++)
	for (w=0; w<ddN.W; w++) 
	  if ( ddS.phi[t][w]<=0 )
	    yap_message(" saving ddS.phi[%d][%d]<=0\n", t, w);
    } 
#endif
  } else if ( ESTIMATE_BETA ) {
    /*
     *  not in memory so have to read
     */
    float *vec;
    FILE *fpin = fopen(phi_file,"rb");
    if ( !fpin )
      yap_sysquit("Cannot open '%s' for read in phi_save()\n",phi_file);
    vec = fvec(ddN.W);
    if ( fseek(fpin, ddN.W*ddN.T*sizeof(vec[0]), SEEK_SET) ) 
      yap_sysquit("Cannot seek in '%s' in phi_save()\n", 
		  phi_file);
    if ( fread(vec, sizeof(vec[0]), ddN.W, fpin) !=ddN.W )
      yap_sysquit("Cannot read vector from '%s' in phi_save()\n", 
		  phi_file);
    fclose(fpin);
    write_fvec(beta_file,ddN.W,vec);
    free(vec);
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
 
void phi_load(char *resstem) {
  uint32_t cnt, dim, size;
  int totgot, got;
  FILE *fp=NULL;
  phi_file = yap_makename(resstem,".phi");
  ddP.phi = fmat(ddN.T+(ESTIMATE_BETA?1:0),ddN.W);
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
  size = ddN.W*(ddN.T+(ESTIMATE_BETA?1:0));
  while ( (got=(int)fread(&ddP.phi[0][got], sizeof(ddS.phi[0][0]), size-totgot, fp)) 
       !=size-totgot ) {
    yap_message("load_phi(%d,%d): got %d in %d\n", (int)ddN.W, (int)ddN.T, got, totgot);
    if ( got==0 || ferror(fp) || feof(fp) )
       yap_sysquit("Cannot read matrix from '%s' in phi_load()\n", phi_file);
    totgot += got;
  }
  totgot += got;
  assert(totgot==size);
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
  for (t=0; t<ddN.T+1; t++) {
    double totvec = 0;
    if ( !ESTIMATE_BETA && t>=ddN.T ) 
      break;
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
      double val = (t>=ddN.T)? betabasewordprob(w) : wordprob(w, t);
      totvec += val;
      if ( val<=0 )
	yap_message("phi_update for phi[%d][%d] zero\n",t,w);
      vec[w] = (ddG.phi_cnt*vec[w] + val) / (ddG.phi_cnt+1);
    }
#ifndef NDEBUG
    if ( t<ddN.T && abs(totvec-1.0)>0.02 )  {
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
  int zerod = 1;
  int t;
  double totvec = 0;
  /*
   *    now compute
   */
  for (t=0; t<ddN.T; t++) {
    double val;
    if ( ddP.PYalpha!=H_HPDD || ddS.TDt[t]>0 || zerod ) {
      val = alphabasetopicprob(t);
      if (zerod) zerod = 0;
    } else
      val = 0;   
    totvec += val;
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
