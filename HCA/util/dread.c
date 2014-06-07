/*
 * Data reading routines
 * Copyright (C) 2010-2013 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *
 *     Data reading only.
 *     They all set:
 *            dbp->N, dbp->W, dbp->D, dbp->d[], dbp->w[]
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>

#include "yap.h"
#include "util.h"
#include "dread.h"

static D_bag_t *dbag_alloc(int N) {
  D_bag_t *dbp = malloc(sizeof(D_bag_t));
  dbp->d = u32vec(N);
  dbp->w = u32vec(N);
  dbp->N = N;
  dbp->D = 0;
  dbp->W = 0;
  return dbp;
}

static D_bag_t *data_read_witdit(char *stem) {
  D_bag_t *dbp;
  int N;
  char *wname, *dname;
  FILE *fp;
  wname = yap_makename(stem, ".wit");
  dname = yap_makename(stem, ".dit");
  /*
   *  file existance test
   */
  fp = fopen(wname,"r");
  if ( !fp )
    yap_sysquit("Cannot open word file '%s'\n", wname);
  fclose(fp);
  fp = fopen(dname,"r");
  if ( !fp )
    yap_sysquit("Cannot open doc file '%s'\n", dname);
  fclose(fp);

  N = countntot(wname);  
  dbp = dbag_alloc(N);
  read_didwid(wname, dname, N, dbp->d, dbp->w);
  dbp->D = u32max(dbp->N,dbp->d); dbp->D++;
  dbp->W = u32max(dbp->N,dbp->w); dbp->W++;
  free(wname);
  free(dname);  
  return dbp;
}

static D_bag_t *data_read_ldac(char *stem) {
  D_bag_t *dbp;
  char *wname;
  FILE  *fp;
  unsigned din, win = 0, nw = 0;
  int i;
  wname = yap_makename(stem, ".ldac");
  fp = fopen(wname,"r");
  if ( !fp )
    yap_sysquit("Cannot open ldac file '%s'\n", wname);
  /*
   *   first run through to get dims
   */
  for (din=0; !feof(fp) && !ferror(fp); din++) {
    unsigned n1, n2, ln;
    if ( fscanf(fp," %u", &ln) != 1 ) {
      if ( feof(fp) )
	break;
      yap_sysquit("Cannot read %d-th data line from '%s'\n", din, wname);
    }
    for (i=0; i<ln; i++) {
      if ( fscanf(fp," %u:%u", &n1, &n2) != 2 )
	yap_sysquit("Cannot read %d-th data line from '%s'\n", ln, wname);
      //  n1--;
      if ( win<n1 )
	win = n1;
      nw += n2;
    }
  }
  if ( ferror(fp) )
    yap_sysquit("error on file '%s'\n", wname);
  dbp = dbag_alloc(nw);
  dbp->D = din;
  dbp->W = win+1;
  yap_message("Read from ldac file: D=%d, W=%d, N=%d\n", dbp->D, dbp->W, dbp->N );
  rewind(fp);
  nw = 0;
  for (din=0; din<dbp->D && !feof(fp); din++) {
    unsigned n1, n2, ln;
    if ( fscanf(fp," %u", &ln) != 1 )
      yap_sysquit("Cannot read %d-th data line from '%s'\n", din, wname);
    for (i=0; i<ln; i++) {
      if ( fscanf(fp," %u:%u", &n1, &n2) != 2 )
	yap_sysquit("Cannot read %d-th entry %d-th data line from '%s'\n", 
		    ln, din, wname);
      // n1--;
      for ( ; n2>0; n2-- ) {
	dbp->d[nw] = din;
	dbp->w[nw++] = n1;
      }
    }
  }
  fclose(fp);
  free(wname); 
  return dbp;
}

static D_bag_t *data_read_docword(char *stem) {
  D_bag_t *dbp;
  char *wname;
  FILE  *fp;
  unsigned din, win, nw, pn, N;
  int i;
  wname = yap_makename(stem, ".docword");
  fp = fopen(wname,"r");
  if ( !fp )
    yap_sysquit("Cannot open docwords file '%s'\n", wname);
  if ( fscanf(fp,"%u %u %u", &din, &win, &nw) != 3 )
    yap_sysquit("Cannot read dimensions from '%s'\n", wname);
  N = 0;
  for (i=0; i<nw; i++) {
    unsigned n1, n2, n3;
    if ( fscanf(fp," %u %u %u", &n1, &n2, &n3) != 3 )
      yap_sysquit("Cannot read %n-th data line from '%s'\n", i, wname);
    N += n3;
    assert(n3>=1);
  }
  assert(N>=nw);
  rewind(fp);
  if ( fscanf(fp,"%u %u %u", &din, &win, &nw)<3 )
    yap_sysquit("Cannot read dimensions from '%s'\n", i, wname);
  dbp = dbag_alloc(N);
  dbp->D = din;
  dbp->W = win;
  for (i=0, pn=0; i<nw; i++) {
    unsigned n1, n2, n3;
    if ( fscanf(fp," %u %u %u", &n1, &n2, &n3) != 3 )
      yap_sysquit("Cannot read %n-th data line from '%s'\n", i, wname);
    for ( ; n3>0 ; n3--) {
      /*
       *  file has indexes offset at 1, convert to offset at 0
       */
      dbp->d[pn] = n1-1;
      dbp->w[pn] = n2-1;
      pn++;
    }
  }
  assert(pn==N);
  fclose(fp);
  free(wname); 
  return dbp;
}

static D_bag_t *data_read_bag(char *stem) {
  D_bag_t *dbp;
  char *wname;
  FILE  *fp;
  unsigned din, win, pn, N;
  int i;
  wname = yap_makename(stem, ".txtbag");
  fp = fopen(wname,"r");
  if ( !fp )
    yap_sysquit("Cannot open txtbag file '%s'\n", wname);
  if ( fscanf(fp,"%u %u", &din, &win) != 2 )
    yap_sysquit("Cannot read dimensions from '%s'\n", wname);
  N = 0;
  for (i=0; i<din; i++) {
    unsigned nl, w, c;
    if ( fscanf(fp," %u", &nl) != 1 )
      yap_sysquit("Cannot read %n-th data line from '%s'\n", i, wname);
    while ( nl-->0 ) {
      if ( fscanf(fp," %u %u", &w, &c) != 2 )
	yap_sysquit("Cannot read %n-th data line from '%s'\n", i, wname);
      N += c;
      if ( w>=win )
        yap_quit("Txtbag file '%s' has word index too large\n", wname);
    }
  }
  rewind(fp);
  if ( fscanf(fp,"%u %u", &din, &win)<2 )
    yap_sysquit("Cannot read dimensions from '%s'\n", i, wname);
  dbp = dbag_alloc(N);
  dbp->D = din;
  dbp->W = win;
  for (i=0, pn=0; i<din; i++) {
    unsigned nl, w, c;
    if ( fscanf(fp," %u", &nl) != 1 )
      yap_sysquit("Cannot read %n-th data line from '%s'\n", i, wname);
    while ( nl-->0 ) {
      if ( fscanf(fp," %u %u", &w, &c) != 2 )
	yap_sysquit("Cannot read %n-th data line from '%s'\n", i, wname);
      while ( c-->0 ) {
	dbp->d[pn] = i;
	dbp->w[pn] = w;
	pn++;
      }
    }
  }
  fclose(fp);
  free(wname); 
  return dbp;
}

static D_bag_t *data_read_lst(char *stem) {
  D_bag_t *dbp;
  int N;
  char *wname;
  FILE  *fp;
  unsigned din, win, pn;
  int i;
  wname = yap_makename(stem, ".txtbag");
  fp = fopen(wname,"r");
  if ( !fp )
    yap_sysquit("Cannot open list txtbag file '%s'\n", wname);
  if ( fscanf(fp,"%u %u", &din, &win) != 2 )
    yap_sysquit("Cannot read dimensions from '%s'\n", wname);
  N = 0;
  for (i=0; i<din; i++) {
    unsigned nl, w;
    if ( fscanf(fp," %u", &nl) != 1 )
      yap_sysquit("Cannot read %n-th data line from '%s'\n", i, wname);
    while ( nl-->0 ) {
      if ( fscanf(fp," %u", &w) != 1 )
	yap_sysquit("Cannot read %n-th data line from '%s'\n", i, wname);
      N ++;
    }
  }
  rewind(fp);
  if ( fscanf(fp,"%u %u", &din, &win)<2 )
    yap_sysquit("Cannot read dimensions from '%s'\n", i, wname);
  dbp = dbag_alloc(N);
  dbp->D = din;
  dbp->W = win;
  for (i=0, pn=0; i<din; i++) {
    unsigned nl, w;
    if ( fscanf(fp," %u", &nl) != 1 )
      yap_sysquit("Cannot read %n-th data line from '%s'\n", i, wname);
    while ( nl-->0 ) {
      if ( fscanf(fp," %u", &w) != 1 )
	yap_sysquit("Cannot read %n-th data line from '%s'\n", i, wname);
      dbp->d[pn] = i;
      dbp->w[pn] = w;
      pn++;
    }
  }
  fclose(fp);
  free(wname); 
  return dbp;
}

char *data_name(char *stem, enum dataType data) {
  char *tname;  
  if ( data==WitDit ) 
    tname = yap_makename(stem,".wit");
  else if ( data==LdaC ) 
    tname = yap_makename(stem,".ldac");
  else if ( data==Docword ) 
    tname = yap_makename(stem,".docword");
  else
    tname = yap_makename(stem,".txtbag");
  return tname;
}

D_bag_t *data_read(char *stem, enum dataType data) {
  D_bag_t *dbp = NULL;
  if ( data==WitDit ) 
    dbp = data_read_witdit(stem);
  else if ( data==LdaC ) 
    dbp = data_read_ldac(stem);
  else if ( data==Docword ) 
    dbp = data_read_docword(stem);
  else if ( data==TxtBag ) 
    dbp = data_read_bag(stem);
  else 
    dbp = data_read_lst(stem);
  return dbp;
}

void data_vocabshrink(D_bag_t *dbp, int maxword) {
  int n, tn;
  if ( dbp->W<=maxword )
    return;
  for (tn=n=0; n<dbp->N; n++) {
    if ( dbp->w[n]<maxword ) {
      dbp->w[tn] = dbp->w[n];
      dbp->d[tn] = dbp->d[n];
      tn++;
    }
  }
  dbp->W = maxword;
  if ( n==tn )
    return;
  if ( tn==0 )
    yap_quit("Shrinking data set to empty!\n");
  dbp->N = tn;
  dbp->d = realloc(dbp->d,sizeof(dbp->d[0])*tn);
  dbp->w = realloc(dbp->w,sizeof(dbp->w[0])*tn);
}

void data_shrink(D_bag_t *dbp, int size) {
  int n;
  if ( size>=dbp->D )
    return;
  for (n=0; n<dbp->N; n++)
    if ( dbp->d[n]>=size ) break;
  if ( n>=dbp->N )
    return;
  if ( n==0 )
    yap_quit("Shrinking data set to empty!\n");
  dbp->N = n;
  dbp->D = size;
  dbp->d = realloc(dbp->d,sizeof(dbp->d[0])*n);
  dbp->w = realloc(dbp->w,sizeof(dbp->w[0])*n);
}

void data_append(D_bag_t *dbp, D_bag_t *dbp2) {
  int i, I, D;
  I = dbp->N;
  D = dbp->D;
  dbp->N += dbp2->N;
  dbp->D += dbp2->D;
  if ( dbp->W < dbp2->W )
    dbp->W = dbp2->W;
  dbp->d = realloc(dbp->d, dbp->N*sizeof(dbp->d[0]));
  dbp->w = realloc(dbp->w, dbp->N*sizeof(dbp->w[0]));
  for (i=I; i<dbp->N; i++) {
    dbp->d[i] = D+dbp2->d[i-I];
    dbp->w[i] = dbp2->w[i-I];
  }
}
