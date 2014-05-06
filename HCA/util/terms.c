/**
 * Main driver
 * Copyright (C) 2011-2013 Wray Buntine
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
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "yap.h"
#include "util.h"
#include "dread.h"

int verbose = 0;
static uint32_t *data_df(char *stem, int *n_df) {
  char buf[1000];
  char *wname = yap_makename(stem, ".srcpar");
  FILE *fp;
  int i;
  int W, D;
  uint32_t *df;
  /*
   *  check .srcpar file exists
   */
  fp = fopen(wname,"r");
  if ( !fp )
    yap_quit("Parameter file '%s' doesn't exist\n", wname);
  fclose(fp);
  free(wname);
  /*  read dfdocs */
  {
    char *p;
    D = atoi(readsrcpar(stem,"documents",buf,50));
    W = atoi(readsrcpar(stem,"features",buf,50));
    if ( D<=0 || W<=0 )
      yap_quit("Cannot read documents or features from .srcpar\n");
    p = readsrcpar(stem,"dfdocs",buf,50);
    if ( p )
      *n_df = atoi(p);
    else
      *n_df = D;
  }
  df = u32vec(W);
  /*  read dfs */
  wname = yap_makename(stem, ".words");
  fp = fopen(wname ,"r");
  if ( !fp )
    yap_sysquit( "Cannot open file '%s' for read\n", wname);
  for (i = 0; i < W; i++) {
    int sl;
    unsigned thisdf;
    if ( fgets(&buf[0],sizeof(buf)-1,fp)==NULL )
      yap_sysquit("Cannot read line %d from '%s'\n", i+1, wname);
    sl = strlen(&buf[0]);
    if ( ! iscntrl(buf[sl-1]) )
      /*   line too long  */
      yap_quit("Cannot parse line %d from '%s', too long\n", i, wname);
    if ( sscanf(&buf[0],"%*u %*s %*x %*u %u ", &thisdf) != 1 )
      yap_quit("Cannot parse line %d from '%s', no df\n", i, wname);
    df[i] = thisdf;
  }
  fclose(fp);
  free(wname);
  return df;
}

/*
 *  copied from Wikipedia page Okapi_BM25
 */
static double bm25(int dl, int tf, int df, double avgdl, int n_df) {
  double k1 = 1.6;
  double b = 0.75;
  double score;
  score = log ((n_df - df + 0.5)/(df + 0.5));
  score *= tf * (k1+1);
  score /= (tf  + k1*(1 - b + b*dl/avgdl));
  return score;
}

void addbest(uint32_t *topw, double *topscore, int *top, int maxtop,
             int w, double score) {
  int ptr2, ptr = *top;
  /*  where to place it */
  while ( ptr>0 && topscore[ptr-1]<score ) {
    ptr--;
  }
  if ( ptr>=maxtop )
    return;

  /*  shuffle */
  ptr2 = *top;
  if ( ptr2>=maxtop )
    ptr2--;
  for ( ; ptr2>ptr; ptr2-- ) {
    topw[ptr2] = topw[ptr2-1];
    topscore[ptr2] = topscore[ptr2-1];
  }
  topw[ptr] = w;
  topscore[ptr] = score;
  if ( *top < maxtop )
    (*top)++;
#if 0 
  {
    int l;
    fprintf(stdout, "addbest:\n");
    for (l=0; l<*top; l++) {
      fprintf(stdout, "%d %lf\n", topw[l], topscore[l]);
    }
  }
#endif
}

/*==========================================
 * main
 *========================================== */
int main(int argc, char* argv[])
{
  enum dataType data = LdaC;
  D_bag_t *dbp;
  char **tokens;
  char *stem;
  int n_df;
  uint32_t *df;
  int i, l;
  int count = 20;
  int MINDF = 5;
  double avdl;
  int *dl;
  int d, w, tf, c;

  uint32_t *topw;
  double *topscore;
  int top=0;

  while ( (c=getopt(argc, argv,"c:f:m:v"))>=0 ) {
    switch ( c ) {
    case 'c':
      if ( !optarg || sscanf(optarg,"%d",&count)!=1 )
	yap_quit("Need a valid 'c' argument\n");
      break;
    case 'f':
      if ( strcmp(optarg,"witdit")==0 ) 
	data = WitDit;
      else if ( strcmp(optarg,"docword")==0 ) 
	data = Docword;
      else if ( strcmp(optarg,"ldac")==0 ) 
	data = LdaC;
      else if ( strcmp(optarg,"bag")==0 ) 
	data = TxtBag;
      else
	yap_quit("Illegal data type for -f\n");
      break;
    case 'm':
      if ( !optarg || sscanf(optarg,"%d",&MINDF)!=1 )
	yap_quit("Need a valid 'm' argument\n");
      break;
    case 'v':
      verbose++;
      break;
    default:
      yap_quit("Unknown option '%c'\n", c);
    }
  }

  if (argc-optind != 1) {
    yap_quit("No file stem\n", c);
    exit(-1);
  }
  stem = strdup(argv[optind++]);
  dbp = data_read(stem,data);

  {
    char *wname = yap_makename(stem, ".tokens");
    tokens = read_vocab(wname,dbp->W,50);
    free(wname);
  }
  
  df = data_df(stem,&n_df);
  dl = u32vec(dbp->D);
  for (i=0; i<dbp->N; i++) {
    dl[dbp->d[i]]++;
  }
  avdl = dbp->N/((double)dbp->D);

  {
    char *dlname = yap_makename(stem,".dl");
    FILE *fp = fopen(dlname,"w");
    for (i=0; i<dbp->D; i++) 
      fprintf(fp,"%d %d\n", i, dl[i]);
    fclose(fp);
    free(dlname);
  }

  topw = u32vec(count);
  topscore = dvec(count);

  w = dbp->w[0];
  d = dbp->d[0];
  tf = 0;
  top = 0;
  for (i=0; i<dbp->N; i++) {
    if ( d == dbp->d[i] ) {
      if ( w == dbp->w[i] ) 
        tf++;
      else if ( tf>0 ) {
        /*  term changed */
        if ( df[w]>MINDF )
          addbest(topw, topscore, &top, count, w, 
                  bm25(dl[d], tf, df[w], avdl, n_df));
        w = dbp->w[i];
        tf = 1;
      } 
    } else {
      /* doc changed */
      if ( tf>0 ) {
        if ( df[w]>MINDF )
          addbest(topw, topscore, &top, count, w, 
                  bm25(dl[d], tf, df[w], avdl, n_df));
        for (l=0; l<top; l++) {
          printf("%d %lf %s\n", d, topscore[l], tokens[topw[l]]);
        }
      }
      top = 0;
      d = dbp->d[i];
      w = dbp->w[i];
      tf = 1;
    }
  }
  if ( tf>0 ) {
    if ( df[w]>MINDF )
      addbest(topw, topscore, &top, count, w, 
              bm25(dl[d], tf, df[w], avdl, n_df));
    for (l=0; l<top; l++) {
      printf("%d %lf %s\n", d, topscore[l], tokens[topw[l]]);
    }
  }

  free(stem);
  free(df);
  return 0;
}
