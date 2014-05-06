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

#define HASHPRIME 7919U
#define REHASHPRIME 7883U

/*
 *   stored as +1 to make 0 be "empty"
 */
static uint32_t *hashtab = NULL;
static uint32_t hashsize = 0;
/*
 *    w = word to find
 *    wind[] = index of words
 */
static int32_t findw(uint32_t w, uint32_t *wind) {
  uint32_t key = ( w*HASHPRIME ) % hashsize;
  if ( hashtab[key]==0 ) 
    return UINT32_MAX;
  while ( wind[hashtab[key]-1] != w ) {
    key = (key+REHASHPRIME) % hashsize;
    if ( hashtab[key]==0 ) 
      return UINT32_MAX;
  }
  return hashtab[key]-1;
}
static void addw(uint32_t w, uint32_t ind) {
  uint32_t key = ( w*HASHPRIME ) % hashsize;
  while ( hashtab[key] != 0 ) {
    key = (key+REHASHPRIME) % hashsize;
  }
  hashtab[key] = ind+1;
  return;
}

/*
 *  print out the topic topk=10 words. report the PMI score. 
 */
double report_pmi(char *topfile,   /* name of topics file */
		  char *pmifile,  /* name of PMI file */
		  int T,          /* total topics */
		  int W,          /* total words */
		  int topk,
		  double *tp)
{
  int lineno = 0;
  int i,k;
  double coherency;
  /*
   *   mapping from local index to actual word index
   */
  uint32_t *wind = u32vec(topk*T);
  int n_wind = 0;
  /*
   *   boolean vector ... is word used
   */
  uint32_t *wuse = u32vec(W/32+1);
  /*
   *   topic words for each topic, using local indices
   */
  uint32_t **topics = u32mat(T,topk);
  /*
   *  PMI's by local index
   */
  double **pmi;

  char *line;
  size_t n_line;
  FILE *fr;
  if ( !wind || !wuse || !topics )
    yap_quit("Out of memory in report_pmi()\n");

  /*
   *   read in file of top word indices in topic
   */
  fr = fopen(topfile,"r");
  if ( !fr ) 
    yap_sysquit("Topic file '%s' not read\n", topfile);
  
  line = NULL;
  n_line = 0;
  while ( getline(&line, &n_line, fr)>0 ) {
    char *buf = line;
    unsigned j;
    lineno ++;
    buf += strspn(buf," \t\n");    //   skip space
    if ( sscanf(buf, "%d: ", &k) <1 )
      yap_quit("Cannot read topic in topic line %d from file '%s'\n", lineno, topfile);
    if ( k<0 || k>=T )
      continue;
    for (i = 0; i<topk && *buf; i++) {
      buf = strpbrk(buf," \t\n");    //   skip to next space
      if ( sscanf(buf, " %u", &j) <1 ) {
	yap_message("Cannot read word %d in topic line %d from file '%s'\n", i+1, lineno, topfile);
	break;
      }
      if ( j>=W) {
	yap_quit("Bad word %d in topic line %d from file '%s'\n", i+1, lineno, topfile);
      }
      buf += strspn(buf," \t\n");    //   skip space
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
	  yap_quit("Lookup of word %d failed at line %d in report_pmi()\n", lineno, (int)j);
	topics[k][i] = ii;
      } else {
	// no, so add it
	wuse[j/32U] |= (1U<<(j%32U));
	wind[n_wind] = j;	
	topics[k][i] = n_wind;
	n_wind++;
      }
    }
    free(line);
    line = NULL;
    n_line = 0;
    if ( i<topk ) topics[k][i] = W;   //  terminator
  }
  fclose(fr);

  pmi = dmat(n_wind,n_wind);
  /*
   *  build hash table now since we know size
   */
  hashsize = n_wind*2;
  hashtab = malloc(sizeof(*hashtab)*hashsize);
  if ( !pmi || !hashtab )
    yap_quit("Out of memory in report_pmi()\n");
  for (i=0; i<hashsize; i++)
    hashtab[i] = 0;
  for (i=0; i<n_wind; i++)
    addw(wind[i],i);

  /*
   *   load up PMI file, only keeping words mentioned in hash table
   */
  {
    unsigned t1, t2;
    double value;
    int zcat = 0;
    fr = fopen(pmifile,"r");
    if ( !fr ) {
      /*
       *    try to zcat it
       */
      char *cmd = malloc(strlen(pmifile)+20);
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
	i1 = findw(t1,wind);
	i2 = findw(t2,wind);
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

  yap_message("PMI: ");
  coherency = 0;
  for (k=0; k<T; k++) {
    double coh = 0;
    int cnt = 0;
    int j;
    for (i=0; i<topk && topics[k][i]<W; i++) {
      for (j=i+1; j<topk && topics[k][j]<W; j++) {
	coh += pmi[topics[k][i]][topics[k][j]];
	cnt ++;
      }
    }
    if ( cnt>0 ) coh /= cnt;
    coherency += coh * tp[k];
    yap_message(" %d:%.3lf", k, coh);
  }
  yap_message(" -> %.3lf\n", coherency);
  
  free(wind);
  free(wuse);
  free(topics[0]); free(topics);
  free(pmi[0]); free(pmi);
  free(hashtab);
  hashtab = NULL;
  hashsize = 0;
  return coherency;
}
