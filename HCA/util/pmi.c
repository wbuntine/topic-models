/*
 * PMI computation from top topics previously computed
 * Copyright (C) 2013 Wray Buntine
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine
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
#include "ehash.h"
extern int verbose;


/****************************************
 *  reading of top topics file
 */
/*
 *    get next (topic,epoch) pair
 *    (T,E) are input maximum topics and epochs
 *    return zero if none
 */
static int lineno;      /* no. lines read */
static char *line;      /*  text of current line */
static char *buf;       /*  ptr to next spot in line */
static char *ttop_file;
static FILE *ttop_fp;

void ttop_eol() {
  if ( line ) {
    free(line);
    line = NULL;
    buf = NULL;
  }
}
void ttop_close() {
  fclose(ttop_fp);
  ttop_file = NULL;
  ttop_fp = NULL;
  ttop_eol();
}      
void ttop_open(char *topfile) {
  ttop_file = topfile;
  ttop_fp = fopen(topfile,"r");
  if ( !ttop_fp ) 
    yap_sysquit("Topic file '%s' not read\n", topfile);
  lineno = 0;
}
int ttop_next(int T, int E, int *k, int *e) {
  size_t n_line = 0;
  line = NULL;
  while ( getline(&line, &n_line, ttop_fp)>0 ) {
    buf = line;
    *e = 0;
    *k = 0;
    lineno ++;
    /*  skip lines starting wth # */
    if ( buf[0]=='#' ) 
      continue;
    buf += strspn(buf," \t\n");    //   skip space
    if ( (E==1 && sscanf(buf, "%d: ", k)<1) || 
	 (E>1 && sscanf(buf, "%d,%d: ", e, k)<2) ) 
      yap_quit("Cannot read topic in topic line %d from file '%s'\n", 
	       lineno, ttop_file);
    /*  skip lines with topics/epochs outside bounds */
    if ( *k<0 || *k>=T )
      continue;
    if ( *e<0 || *e>=E )
      continue;
    return 1;
  } 
  if ( line ) {
    free(line);
    line = NULL;
  }
  return 0;
}
int ttop_word(int W, int n_word, unsigned int *j) {
  while ( *buf ) {
    buf = strpbrk(buf," \t\n");    //   skip to next space
    if ( sscanf(buf, " %u", j) <1 ) {
      if ( verbose>2 ) 
        yap_message("Cannot read word %d in topic line %d from file '%s'\n", 
		    n_word+1, lineno, ttop_file);
      break;
    }
    if ( *j>=W) {
      yap_quit("Bad word %d in topic line %d from file '%s'\n", 
               n_word+1, lineno, ttop_file);
    }
    buf += strspn(buf," \t\n");    //   skip space
    return 1;
  }
  free(line);
  line = NULL;
  buf = NULL;
  return 0;
}
/*
 *    if tpmi==NULL
 *         print out PMI for topics computed on topk words
 *    else
 *         store in tpmi
 *    tpmi[e*(T+1)+k] = PMI for topic k in epoch e
 *    tpmi[e*(T+1)+T] = mean (by topic propportion) PMI for epoch 
 *
 *    return average over epochs
 */
double report_pmi(char *topfile,  /* name of topics file to read */
		  char *pmifile,  /* name of PMI file to read */
                  char *toppmifile, /* name of topics+pmi file to write */
		  int T,          /* total topics */
		  int W,          /* total words */
		  int E,          /* number of epochs */
		  int topk,
		  double *tp,
                  float *tpmi)
{
  int k, e, thee;
  /*
   *   mapping from local index to actual word index
   */
  uint32_t *wind = u32vec(topk*T*E);
  int n_wind = 0;
  /*
   *   boolean vector ... is word used
   */
  uint32_t *wuse = u32vec(W/32+1);
  /*
   *  PMI's by local index
   */
  uint32_t *topic = u32vec(topk);
  float *coherency = fvec(E);
  double **pmi;
  float ave = 0;
  FILE *tpfp = NULL;
  ehash_t hp;

  if ( !wind || !wuse )
    yap_quit("Out of memory in report_pmi()\n");

  /*
   *   read in file of top word indices in topic, but
   *   all we do first is collect all the words
   */
  ttop_open(topfile);
  while ( ttop_next(T, E, &k, &e) ) {
    int i;
    unsigned j;
    for (i = 0; i<topk && ttop_word(W, i, &j); i++) {
      /*
       *   check if word exists, and set up its index
       */
      if ( wuse[j/32U] & (1U<<(j%32U)) ) {
	// yes, so search for it
	int ii;
	for (ii=0; ii<n_wind; ii++)
	  if ( wind[ii]==j )
	    break;
	if ( ii>=n_wind )
	  yap_quit("Lookup of word %d failed at line %d in report_pmi()\n", 
		   (int)j, lineno);
      } else {
	// no, so add it
	wuse[j/32U] |= (1U<<(j%32U));
	wind[n_wind] = j;	
	n_wind++;
      }
    }
    ttop_eol();
  }
  ttop_close();

  pmi = dmat(n_wind,n_wind);
  if ( !pmi )
    yap_quit("Out of memory in report_pmi()\n");
  /*
   *  build hash table now since we know size
   */
  if ( ehash_init(&hp, n_wind*2) )
    yap_quit("Out of memory in report_pmi()\n");
  for (k=0; k<n_wind; k++)
    ehash_addw(&hp,wind[k],k);

  /*
   *   load up PMI file, only keeping words mentioned in hash table
   */
  {
    unsigned t1, t2;
    double value;
    int zcat = 0;
    FILE *fr;
    fr = fopen(pmifile,"r");
    if ( !fr ) {
      /*
       *    try to zcat it
       */
      char *cmd = malloc(strlen(pmifile)+20);
      sprintf(cmd,"%s.gz", pmifile);
      fr = fopen(cmd,"r");
      if ( !fr ) 
	yap_sysquit("Cannot open pmifile '%s' in report_pmi()\n", 
		    pmifile);
      fclose(fr);
      sprintf(cmd,"gunzip -c %s", pmifile);
      fr = popen(cmd,"r");
      if ( !fr )
	yap_sysquit("Cannot open or zcat pmifile '%s' in report_pmi()\n", 
		    pmifile);
      zcat = 1;
      free(cmd);
    }
    while (fscanf(fr, "%u %u %lg", &t1, &t2, &value)==3 ) { 
      if ( t1>=W || t2>= W )
	yap_quit("Illegal word index in report_pmi()\n");
      if ( t1!= t2 && ( wuse[t1/32U] & (1U<<(t1%32U)) ) 
	   && ( wuse[t2/32U] & (1U<<(t2%32U))) ) {
	int i1, i2;
	i1 = ehash_findw(&hp, t1, wind);
	i2 = ehash_findw(&hp, t2, wind);
	if ( i1==UINT32_MAX || i2==UINT32_MAX )
	  yap_quit("Could not locate word index in report_pmi()\n");
	pmi[i1][i2]=value;
	pmi[i2][i1]=value;
      }
    }
    if ( zcat )
      pclose(fr);
    else
      fclose(fr);
  }
  
  /*
   *    compute PMI score for each topic
   */
  ttop_open(topfile);
  thee = 0;
  if ( tpmi==NULL ) {
    if ( E>1 ) 
      yap_message("PMI %d:: ", 0);
    else
      yap_message("PMI :: ");
  }
  if ( toppmifile ) {
    tpfp = fopen(toppmifile,"w");
    if ( !tpfp ) 
      yap_sysquit("Cannot open '%s' for write\n", toppmifile);
  }
  while ( ttop_next(T, E, &k, &e) ) {
    int i;
    unsigned j;
    int cnt = 0;
    int e = 0;
    double coh = 0;
    if ( e!=thee ) {
      thee = e;
      if ( tpmi==NULL )
        yap_message("\nPMI %d:: ", e);
    }
    for (i = 0; i<topk && ttop_word(W, i, &j); i++) {
      topic[i] = ehash_findw(&hp, j,wind);
    }
    ttop_eol();
    if ( i<topk )
      topic[i] = W;
    /*
     *  topics now read 
     */
    if ( tpfp ) {
      if ( E>1 )
        fprintf(tpfp, "%d,%d: ", e, k);
      else
        fprintf(tpfp, "%d: ", k);
    }
    for (i=0; i<topk && topic[i]<W; i++) {
      double cohone = 0;
      int cntone = 0;
      for (j=0; j<topk && topic[j]<W; j++) {
	cohone += pmi[topic[i]][topic[j]];
	cntone ++;
      }
      coh += cohone;
      cnt += cntone;
      if ( tpfp ) 
        fprintf(tpfp, " %d(%lf)", wind[topic[i]], cohone/cntone);
    }
    coh /= 2;
    if ( cnt>0 ) coh /= cnt;
    coherency[e] += coh * tp[k];
    if ( tpmi )
      tpmi[e*(T+1)+k] = coh;
    else
      yap_message(" %d:%.3lf", k, coh);
    if ( tpfp )
      fprintf(tpfp, " -> %lf(%lf)\n", coh, tp[k]);
  }
  ttop_close();
  if ( tpfp )
    fclose(tpfp);

  if ( tpmi==NULL ) yap_message("\nPMI =");
  if ( E==1 ) {
    if ( tpmi==NULL ) yap_message(" %.3lf\n", coherency[0]);
    ave = coherency[0];
    tpmi[T] = ave;
  } else {
    int e;
    for (e=0; e<E; e++) {
      ave += coherency[e];
      if ( tpmi ) 
        tpmi[e*(T+1)+T] = coherency[e];
      else
        yap_message(" %.3lf", coherency[e]);
    }
    ave /= E;
    if ( tpmi==NULL ) yap_message(" -> %.3lf\n", ave);
  }
      
  free(wind);
  free(coherency);
  free(wuse);
  free(topic);
  free(pmi[0]); free(pmi);
  ehash_free(&hp);
  return ave;
}
