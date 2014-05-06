/*
 * Module for bursty data preparation and processing
 * Copyright (C) 2013-4 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *
 *   dmi_*() routines prepare the collection:
 *           need to record which words occur multiple times in docs
 *   misi_*() routines do per document processing
 *           see dmi_rand() for usage
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
#include "srng.h"
#include "misi.h"
#include "gibbs.h"

// #define MISI_CHECK

static int zero(int l) {
  return 0;
}

/*
 *  allocation
 */
void misi_init(D_DMi_t *pdmi, D_MiSi_t *ptr) {
  ptr->pdmi = pdmi;
  ptr->Mik = u16mat(pdmi->MI_max, pdmi->T);
  ptr->Sik = u16mat(pdmi->MI_max, pdmi->T);
  ptr->Mi = u16vec(pdmi->T);
  ptr->Si = u16vec(pdmi->T);
}

void misi_free(D_MiSi_t *ptr) {
    free(ptr->Mik[0]); free(ptr->Mik); free(ptr->Mi);
    free(ptr->Sik[0]); free(ptr->Sik); free(ptr->Si);
}

/*
 *   call to ensure stats are zero for a document, for
 *   use in LRS
 */
void misi_zero(D_MiSi_t *ptr, int d) {
  int t, mi;
  D_DMi_t *pdmi = ptr->pdmi;
  for (t=0; t<pdmi->T; t++)
    ptr->Mi[t] = ptr->Si[t] = 0;
  /*  figure out lowest pdmi->multiind[mi] for this doc */
  mi = pdmi->MI[d];
  ptr->mi_base = pdmi->multiind[mi];
  for ( ; mi<pdmi->MI[d+1]; mi++)
    if ( ptr->mi_base>pdmi->multiind[mi] )
      ptr->mi_base = pdmi->multiind[mi];
}

/*
 *   build all stats for a one document
 */
void misi_build(D_MiSi_t *ptr, int d, int notest) {
  int t, l, mi;
  D_DMi_t *pdmi = ptr->pdmi;
  misi_zero(ptr, d);
  /*  compute Mi[] and Si[] */
#ifndef NDEBUG 
  if ( !notest ) {
    for (mi=0; mi<pdmi->MI_max; mi++) {
      for (t=0; t<pdmi->T; t++) {
	assert(ptr->Sik[mi][t]==0);
	assert(ptr->Mik[mi][t]==0);
      }
    }   
  }
#endif
  mi = pdmi->MI[d];      
  for (l=pdmi->NdTcum[d]; l< pdmi->NdTcum[d+1]; l++) {
    if ( !pdmi->holdout(l) ){
      t = Z_t(pdmi->z[l]);
      ptr->Mi[t]++;
      if ( misi_multi(pdmi,l) ) {
	int mii = pdmi->multiind[mi] - ptr->mi_base;
	ptr->Mik[mii][t]++;
	if ( Z_issetr(pdmi->z[l]) ) {
	  ptr->Sik[mii][t]++;
	  ptr->Si[t]++;
	}
      } else 
	ptr->Si[t]++;
    }
    if ( misi_multi(pdmi,l) ) 
      mi++;
  }
#ifndef NDEBUG 
  if ( !notest ) {
    for (mi=0; mi<pdmi->MI_max; mi++) {
      for (t=0; t<pdmi->T; t++) {
	assert(ptr->Sik[mi][t]<=ptr->Mik[mi][t]);
	assert(ptr->Mik[mi][t]==0 || ptr->Sik[mi][t]>0);
      }
    }   
  }
#endif
}

/*
 *   shares logic with misi_build() to zero arrays
 */
void misi_unbuild(D_MiSi_t *ptr, int d, int notest) {
  int t, l, mi;
  D_DMi_t *pdmi = ptr->pdmi;
  mi = pdmi->MI[d];      
  for (l=pdmi->NdTcum[d]; l< pdmi->NdTcum[d+1]; l++) {
    if ( misi_multi(pdmi,l) ) {
      if ( !pdmi->holdout(l) ) {
	int mii = pdmi->multiind[mi] - ptr->mi_base;
	t = Z_t(pdmi->z[l]);
	ptr->Mik[mii][t] = 0;
	ptr->Sik[mii][t] = 0;
      }
      mi++;
    }
  }
#ifndef NDEBUG 
  if ( !notest ) {
    for (mi=0; mi<pdmi->MI_max; mi++) {
      for (t=0; t<pdmi->T; t++) {
	assert(ptr->Sik[mi][t]==0);
	assert(ptr->Mik[mi][t]==0);
      }
    }   
  }
#endif
}

int misi_blocked(D_MiSi_t *ptr, int i, int mi, int t) {
  D_DMi_t *pdmi = ptr->pdmi;
  if ( misi_multi(pdmi,i) && Z_issetr(pdmi->z[i]) ) {
    int mii = pdmi->multiind[mi] - ptr->mi_base;
    if ( ptr->Mik[mii][t]>1 && ptr->Sik[mii][t]==1 )
      return 1;
  }
  return 0;
}


void misi_decr(D_MiSi_t *ptr, int i, int mi, int t, int w) {
  D_DMi_t *pdmi = ptr->pdmi;
  if ( misi_multi(pdmi,i)) {
    int mii = pdmi->multiind[mi] - ptr->mi_base;
    assert(ptr->Mik[mii][t]>0);
    if ( w>=0 ) {
      ptr->Sik[mii][t]--; 
    }
    ptr->Mik[mii][t]--; 
    assert(ptr->Mik[mii][t]>=ptr->Sik[mii][t]);
  }
  if ( w>=0 ) {
    ptr->Si[t]--; 
  }
  ptr->Mi[t]--; 
}

int misi_incr(D_MiSi_t *ptr, int i, int mi, int t, int w, float dtip) {
  D_DMi_t *pdmi = ptr->pdmi;
  if ( misi_multi(pdmi,i) ) {
    int mii = pdmi->multiind[mi]-ptr->mi_base;
    if ( dtip==1 || dtip > rng_unit(rngp) )
      Z_setr(pdmi->z[i]);
    else {
      /*   subsequently use -ve wid to flag no word stats */
      w = -1;
      Z_unsetr(pdmi->z[i]);
    }
    if ( w>=0 ) 
      ptr->Sik[mii][t]++; 
    ptr->Mik[mii][t]++; 
  } 
  if ( w>=0 ) 
    ptr->Si[t]++; 
  ptr->Mi[t]++;
  return w;
}


/*
 *   set the first occurrence of a multi to have r=1
 *   but ignore testing needs for now
 */
void dmi_rand(D_DMi_t *pdmi, int firstdoc, int lastdoc) {
  int mi, d, i;
  D_MiSi_t dD;
  misi_init(pdmi,&dD);
  mi = pdmi->MI[firstdoc];
  for (d=firstdoc; d<lastdoc; d++) {
    misi_build(&dD,d,1);  
    for (i=pdmi->NdTcum[d]; i<pdmi->NdTcum[d+1]; i++) {
      if ( !pdmi->holdout(i) ) {
        if ( misi_multi(pdmi,i) ) {
          int mii = pdmi->multiind[mi]-dD.mi_base;
          int t = Z_t(pdmi->z[i]);
          assert(dD.Mik[mii][t]>0);
          if ( dD.Sik[mii][t]==0 && ! Z_issetr(pdmi->z[i]) ) {
            Z_setr(pdmi->z[i]);
            dD.Sik[mii][t]++;
            dD.Si[t]++;
          }
        } else {
          Z_setr(pdmi->z[i]);
        }
      }
      if ( misi_multi(pdmi,i) ) mi++;
    }
    misi_unbuild(&dD,d,0);
  }
  misi_free(&dD);
}

void dmi_check(D_DMi_t *pdmi, int i) {
  D_MiSi_t dD;
  int mi, t;
  int imax = i+1;
  misi_init(pdmi, &dD);
  if ( i<0 ) {
    imax = pdmi->DT;
    i = 0;
  }
  for (; i<imax; i++) {    
    misi_build(&dD,i,1);
    for (mi=0; mi<pdmi->MI_max; mi++) {
      for (t=0; t<pdmi->T; t++) {
	assert(dD.Sik[mi][t]<=dD.Mik[mi][t]);
	assert(dD.Mik[mi][t]==0 || dD.Sik[mi][t]>0);
      }
    }   
    misi_unbuild(&dD,i,1);
  }
  misi_free(&dD);
}

/*
 *
 */
void dmi_init(D_DMi_t *ptr, 
              uint16_t *z, uint32_t *pw, uint32_t *NdTcum, //  cumsum(NdT)
              int T, int N, int W, int D, int DT,   // dims
	      int (*holdout)(int l)) {
  int i, w;
  uint16_t *wcount = (uint16_t*)malloc(sizeof(wcount[0])*W);
  int32_t *wind = (int32_t*)malloc(sizeof(wind[0])*W);
  int dim_Mi = 0;
  int dim_multiind = 0;
  int maxwc = 0;
  if ( !wcount || !wind ) 
    yap_quit("Cannot allocate memory for vecs/matrices\n");
  if ( holdout!=NULL )
    ptr->holdout = holdout;
  else 
    ptr->holdout = &zero;
  ptr->DT = DT;
  ptr->T = T;
  ptr->z = z;
  ptr->NdTcum = NdTcum;
  // a bit vector
  ptr->multi = u32vec((N+31)/32);
  for (w=0; w<W; w++) {
    wind[w] = -1;
    wcount[w] = 0;
  }
  for (i=0; i<((N+31)/32); i++)
    ptr->multi[i] = 0;
  /*
   *   first get dimensions of multiind[] and Mi[];
   *   so go through each document, counting words in wcount[]
   */
  for (i=0; i<D; i++) {
    unsigned l;
    for (l=NdTcum[i]; l<NdTcum[i+1]; l++) {
      wcount[pw[l]]++;
      if ( wcount[pw[l]]==2 ) {
        /*  this word occurs more than once   */
        dim_Mi ++;
        dim_multiind += 2;
      } else if ( wcount[pw[l]]>2 ) {
        dim_multiind ++;
        if ( wcount[pw[l]]>maxwc )
          maxwc = wcount[pw[l]];
      }
    }
    //   set multi[]
    for (l=NdTcum[i]; l<NdTcum[i+1]; l++) 
      if ( wcount[pw[l]]>1 )
        ptr->multi[l/32] |= (1U << (l%32U));
    //   reset wcount[] to zero
    for (l=NdTcum[i]; l<NdTcum[i+1]; l++) 
      wcount[pw[l]] = 0;
    if ( i==DT-1 )
      ptr->dim_MiT = dim_Mi;
  }
  /*
   *   now make and fill out multiind[] and MI
   */
  ptr->dim_multiind = dim_multiind;
  ptr->multiind = u32vec(dim_multiind);
  ptr->MI = u32vec(D+1);
  ptr->Mi = u16vec(dim_Mi);
  ptr->dim_Mi = dim_Mi;
  ptr->Mi_max = maxwc;
  dim_Mi = 0;
  dim_multiind = 0;
  for (i=0; i<D; i++) {
    unsigned l;
    /*  first run through to fill wind[] and .MI[] */
    for (l=NdTcum[i]; l<NdTcum[i+1]; l++) {
      wcount[pw[l]]++;
      if ( wcount[pw[l]]==2 ) {
        ptr->MI[i] += 2;
        wind[pw[l]] = dim_Mi;
        dim_Mi++;
      } else if ( wcount[pw[l]]>2 ) {
        ptr->MI[i] ++;
      }
    }
    /*  now fill ptr->multiind[]  */
    for (l=NdTcum[i]; l<NdTcum[i+1]; l++) {
      if ( wind[pw[l]]>=0 ) 
        ptr->multiind[dim_multiind++] =  wind[pw[l]];
    }
    for (l=NdTcum[i]; l<NdTcum[i+1]; l++) {
      wind[pw[l]] = -1;
      wcount[pw[l]] = 0;
    }
  }
  /*
   *   make MI_max be larger than range of ptr->multiind[] for a doc
   *   NB. computing max_i (ptr->MI[i]) gives range of domain
   *       for ptr->multiind[] , halving it must be upper bound
   */
  ptr->MI_max = 0;
  for (i=0; i<D; i++) 
    if (  ptr->MI_max<ptr->MI[i] )
      ptr->MI_max = ptr->MI[i];
  ptr->MI_max = (ptr->MI_max+1)/2;
  /*
   *   make MI[] index start of ptr->multiind for this doc
   */
  for (i=1; i<D; i++) 
    ptr->MI[i] += ptr->MI[i-1];
  for (i=D; i>0; i--) 
    ptr->MI[i] = ptr->MI[i-1];
  ptr->MI[0] = 0;
  /*
   *   now fill out Mi[] since we have to deal with the holdout
   *   case in a different way:
   *       i.e.,   ptr->multiind[] and ptr->MI[] set up the same
   *               but ptr->Mi[] only gives stats for training 
   *               part of the document
   */
  for (i=0; i<D; i++) {
    unsigned l;
    dim_multiind = ptr->MI[i];
    for (l=NdTcum[i]; l<NdTcum[i+1]; l++) {
      if ( misi_multi(ptr,l) ) {
        if ( !ptr->holdout(l) )
          ptr->Mi[ptr->multiind[dim_multiind]]++;
        dim_multiind++;
      }
    }
  } 
#ifndef NDEBUG
  if ( 0 ) {
    int Mitot = 0;
    int Mimax = 0;
    /*   sanity check  */
    for (i=0; i<ptr->dim_Mi; i++) Mitot += ptr->Mi[i];
    assert(Mitot = ptr->MI[D]);
    for (i=0; i<ptr->dim_Mi; i++)
      if ( ptr->Mi[i]>Mimax )
        Mimax = ptr->Mi[i];
    yap_message("Max of Mi = %d vs %d\n", Mimax, maxwc);
    Mimax = 0;
    for (i=0; i<ptr->dim_multiind; i++)
      if ( ptr->multiind[i]>Mimax )
        Mimax = ptr->multiind[i];
    yap_message("Dim of Mi = %d versus %d\n",
                Mimax+1, ptr->dim_Mi);
  }
#endif
  free(wcount);
  free(wind);
}

void dmi_free(D_DMi_t *ptr) {
  free(ptr->multi);
  free(ptr->multiind);
  free(ptr->Mi);
  free(ptr->MI);
}
