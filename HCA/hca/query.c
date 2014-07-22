/*
 * Query support
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
#include "pctl.h"
#include "probs.h"

/*
 *   a query is a mapping from the word indices to
 *   the position in the query; non-query words map to -1;
 *   the file has multiple lines in format:
 *        NW, W1, W2, ...
 *   where NW = #words, Wk = 0-offset index of word;
 *   so each word assumed to exist only once and ignored otherwise;
 *   each line is one query
 */
#define QMAX 1000
void query_read(char *qname) {
  FILE  *fp;
  unsigned win = 0, nw = 0, nq=0, qin;
  int i;
  uint32_t *wlist = malloc(sizeof(wlist[0])*QMAX);
  uint32_t *qlist = malloc(sizeof(qlist[0])*QMAX);
  int16_t *map = malloc(sizeof(map[0])*ddN.W);
  if ( !map || !wlist || !qlist)
    yap_quit("Cannot allocate memory in query_read()\n");
  fp = fopen(qname,"r");
  if ( !fp )
    yap_sysquit("Cannot open query bag file '%s'\n", qname);
  for (i=0; i<ddN.W; i++) 
    map[i] = -1;
  nw = 0;  qin = 0;
  while ( fscanf(fp," %u", &nw) == 1 ) {
    for (i=0; i<nw; i++) {
      if ( fscanf(fp," %u", &win) != 1 || win>=ddN.W )
	yap_sysquit("Cannot read %d-th entry from '%s'\n", 
		    i, qname);
      if ( map[win]<0 ) {
	qlist[nq] = qin;
	wlist[nq] = win;
	map[win] = nq++;
	if ( nq>=QMAX ) 
	  yap_quit("Predefined query length maximum (%d) too small\n", QMAX);
      } else {
	/*  
	 *    word appears already:  same query, drop, other query, copy
	 */
	if ( qlist[map[win]]!=qin ) {
	  qlist[nq] = qin;
	  wlist[nq] = win;
	  nq++;
	  if ( nq>=QMAX ) 
	    yap_quit("Predefined query length maximum (%d) too small\n", QMAX);
	}
      }
    }
    qin++;
    nw = 0;
  }
  if ( ferror(fp) )
    yap_sysquit("Cannot read data line from '%s'\n", qname);
  fclose(fp);
  ddP.query = map;
  ddP.qid = realloc(qlist, nq*sizeof(qlist[0]));
  ddP.qword = realloc(wlist, nq*sizeof(wlist[0]));
  ddP.qposn = malloc(sizeof(ddP.qposn[0])*(qin+1));
  if ( !ddP.qword || !ddP.qid || !ddP.qposn )
    yap_quit("Cannot allocate memory in query_read()\n");
  ddP.n_words = nq;
  ddP.n_query = qin;
  ddP.qposn[0] = 0;
  for (i=1; i<ddP.n_words; i++) {
    if ( ddP.qid[i] != ddP.qid[i-1] )
      ddP.qposn[ddP.qid[i]] = i;
  }
  ddP.qposn[ddP.n_query] = ddP.n_words;
}

static int besttopic(int w, double *tp) {
  int t, bestt = 0;
  double bestval = wordprob(w,0) * tp[0];
  for (t=1; t<ddN.T; t++) {
    double newval = wordprob(w,t) * tp[t];
    if ( newval>bestval ) {
      bestt = t;
      bestval = newval;
    }
  }
  return bestt;
}

/*
 *   return index where the new entry is placed
 */
static int bubble(int K, int *topind, float *score, float newscore) {
  int k;
  int newind = topind[K-1];
  if ( K>1 ) {
    int scale;
    /*  initial bisection search */
    for (scale=2; scale<K; scale*=2) ;
    scale /= 2;
    k = scale;
    assert(k<K);
    while ( scale>=1  ) {
      assert(topind[k]<K);
      if ( score[topind[k]] > newscore ) {
	k -= scale;
        if ( k<0 ) k = 0;
      } else if ( score[topind[k]] < newscore ) {
	k += scale;
        if ( k>=K ) k = K-1;
      } else
	break;
      scale /= 2;
    }
    /*  make sure not out of bounds */
    if ( k<0 )
      k = 0;
    if ( k>=K )
      k = K-1;
    /*  now do classic bubble, since don't know where landed */
    while ( k<K-1 && score[topind[k]] < newscore ) {
      k++;
    } 
    while ( k>0 && score[topind[k-1]] > newscore ) {
      k--;
    }
    /*
     *  want to finish with score[topind[k-1]] <= newscore and 
     *                      score[topind[k]] >= newscore
     */
    assert(k==0 || score[topind[k-1]] <= newscore );
    assert(score[topind[k]] >= newscore );
    /*
     *    it goes here, ripple down and place
     */
    if ( k<K-1) { 
      int i;
      for (i=K-1; i>k; i--) {
	topind[i] = topind[i-1];
      }
    }
  } else
    k = 0;
  topind[k] = newind;
#if 0 
  score[newind] = newscore;
  for (k=K-1; k>0; k--) {
    assert(score[topind[k-1]]<=score[topind[k]]);
  }
#endif
  return newind;
}

static int n_df;
static uint32_t *df = NULL;

/*
 *  copied from Wikipedia page Okapi_BM25
 */
static double bm25(int d, int *found, uint32_t *wi, int nw,  float *ws) {
  double k1 = 1.6;
  double b = 0.75;
  double rank = 0;
  double avgdl = ((double)ddN.NT)/((double)ddN.DT);
  int j;
  assert(nw>0);
  assert(d>=0 && d<ddN.DT);
  for (j=0; j<nw; j++) {
    double score = log ((n_df - df[wi[j]] + 0.5)/(df[wi[j]] + 0.5));
    score *= found[j] * (k1+1);
    score /= (found[j]  + k1*(1 - b + b*ddD.NdT[d]/avgdl));
    ws[j] = score;
    rank += score;
  }
  return rank;
}

/*
 *   logic taken from core of Gibbs samples
 */
static void query_docprob(int did, int *mimap, float *p, D_MiSi_t *dD,
			  float *cnt, float *wordscore) {
  int l, t, wid;
  double Z, tot;
  int Td_ = 0;
  double *tp;
  
  /*
   *   doing estimation, not sampling so use *prob() versions
   *   of estimates, not *fact() versions
   */

  tp = dvec(ddN.T);
  if ( ddP.PYalpha )
    Td_ = comp_Td(did);
  for (t=0; t<ddN.T; t++) 
    tp[t] = topicprob(did,t,Td_);

  for (l=0; l<ddP.n_words; l++) {
    int cmax = 0;
    wid = ddP.qword[l];
    if ( ddP.query[wid]!=l )
      /*  word has occurred before so drop */
      continue;
    for (t=0, Z=0, tot=0; t<ddN.T; t++) {
      /*
       *   doing estimation, not sampling so use prob versions
       */
      double tf = tp[t];
      if ( tf>0 ) {
	double wf = wordprob(wid, t);
	tot += tf;
	if ( ddP.bdk!=NULL ) {
	  int n, s;
          /*
           *  with burstiness;
	   *  reproduce some logic in docprob() but
	   *  we've got local data structures
           */
	  if ( mimap[l]>ddN.N ) {
	    /*
	     *   doesn't occur in doc
	     */
	    n = s = 0;
	  } else if ( mimap[l]<0 ) {
	    /*
	     *   occurs once in doc
	     */
	    int z = Z_t(ddS.z[-mimap[l]-1]);
	    n = s = (z==t)?1:0;
	  } else {
	    /*
	     *   its a multi
	     */
	    int mii = ddM.multiind[mimap[l]]-dD->mi_base;
	    assert(mii>=0);
	    assert(mii<ddM.MI_max);
	    n = dD->Mik[mii][t];
	    s = dD->Sik[mii][t];
	  } 
	  wf = (wf*(ddP.bdk[t]+ddP.ad*dD->Si[t]) + (n-ddP.ad*s))/
	    (ddP.bdk[t]+dD->Mi[t]); 
	  if ( cmax<n )
	    cmax = n;
	}
	Z += p[t] = tf*wf;
      } else
	p[t] = 0;
    }
    if ( ddP.bdk!=NULL )
      cnt[l] += cmax;
      wordscore[l] += -log(Z/tot);
  }
  free(tp);
}

static void map_query(int d, int *map, int *found) {
  int mi;
  int j, l;
  /*  by default the query has no words in doc d */
  for (j=0; j<ddP.n_words; j++) {
    found[j] = 0;
  }
  if ( map ) {
    mi = ddM.MI[d];
    for (j=0; j<ddP.n_words; j++) {
      map[j] = ddN.N+1;
    }
  }
  for (l=ddD.NdTcum[d]; l<ddD.NdTcum[d+1]; l++) {
    j = ddP.query[ddD.w[l]];
    if ( j>=0 ) {
      found[j]++;
      /*  this word is in the query */
      if ( map ) {
	if ( M_multi(l) )
	  /*   give the index into Mi[] and Si[]  */
	  map[j] = mi;
	else
	  /*  
	   *  tag it as a word occurring once only, and tell
	   *   where the word is
	   */
	  map[j] = -(l+1);
      }
    }
    if ( map && M_multi(l) ) 
      mi++;
  }
}

/*
 *    run regular gibbs cycles on the data with phi used;
 *    the evaluation on each doc, and sample word probs
 *
 *    if qparts>0, split collection into parts and only search this
 *
 *    K = number of top results to retain
 */
void gibbs_query(char *stem, int K, char *qname, int dots, int this_qpart, int qparts) {
  /*
   *    mapping from query word posn. to its mi in current doc
   *       >ddN.N  = not in current doc
   *       -ve  = has no mi since occurs just once, found at
   *              posn  (-map[]-1)
   *       non -ve = mi value
   */
  int     *mimap = NULL;
  /*
   *     usual stuff for Gibbs loop over docs
   */
  int i, j;
  float *fact = fvec(ddN.T*4);
  D_MiSi_t dD;
  /*
   *   an index into topk[] which maintains ordering
   */
  int     *topind;

  /*
   *    these store statistics of the results, for printing
   *    these are unordered, ordered by topind[]
   */
  /*      document score  */
  float   *topscore;
  /*      document number   */
  int     *topk;
  /*      flags if ord is irrelevant, thus not scored  */
  char    *wordunused;

  /*
   *    per word stats for top results saved
   */
  int     *found;
  float   *topcnt;
  float   *topwordscore;
  /*
   *    temporary versions for when gibbs running
   */
  int     *found_buf;
  float   *topcnt_buf;
  float   *topwordscore_buf;
  double  *logprob;
  /*
   *   search here
   */
  int startdoc = 0;
  int enddoc = ddN.DT;

  /*
   *    setup
   */
  topcnt = malloc(sizeof(topcnt[0])*K*ddP.n_words);
  topwordscore = malloc(sizeof(topwordscore[0])*K*ddP.n_words);
  found = malloc(sizeof(found)*ddP.n_words*K);
  wordunused = malloc(sizeof(wordunused[0])*ddP.n_words);

  topcnt_buf = malloc(sizeof(topcnt[0])*ddP.n_words);
  topwordscore_buf = malloc(sizeof(topwordscore[0])*ddP.n_words);
  found_buf = malloc(sizeof(found)*ddP.n_words);
  if ( !topcnt || !topwordscore || !found || 
       !topcnt_buf || !topwordscore_buf || !found_buf )
    yap_quit("Cannot allocate memory in gibbs_query()\n");

  logprob = malloc(sizeof(logprob[0])*ddP.n_query);
  topscore = malloc(sizeof(topscore[0])*K*ddP.n_query);
  topind = malloc(sizeof(topind[0])*K*ddP.n_query);
  topk = malloc(sizeof(topk[0])*K*ddP.n_query);
  if ( ddP.bdk!=NULL ) 
    mimap = malloc(sizeof(mimap[0])*ddP.n_words);
  if ( !topk || !topscore || !logprob || !topind )
    yap_quit("Cannot allocate memory in gibbs_query()\n");
  for (i=0; i<ddP.n_words; i++) {
    wordunused[i] = 0;
  }
  for (i=0; i<K*ddP.n_query; i++) {
    topind[i] = i%K;
    topk[i] = -1;
    topscore[i] = HUGE_VAL;
  }
  
  /*
   *  check words to exclude using topics
   */
  if ( ddP.n_excludetopic>0 ) {
    double *tprob = malloc(sizeof(tprob[0])*ddN.T);
    assert(ddS.Ndt);
    get_probs(tprob);
    yap_probs();
    if ( verbose>1 )
      yap_message("Excluding words: ");
    for (i=0; i<ddP.n_words; i++) {
      int t = besttopic(ddP.qword[i],tprob);
      if ( Q_excludetopic(t) ) {
	wordunused[i] = 1;
	if ( verbose>1 )
	  yap_message(" %d/%d", (int)ddP.qword[i], t);
      }
    } 
    if ( verbose>1 )
      yap_message("\n");
    free(tprob);
  }
  
  if ( ddP.bdk!=NULL ) misi_init(&ddM,&dD);

  if ( qparts>0 ) {
    startdoc = ((double)this_qpart)/qparts * ddN.DT;
    enddoc = ((double)this_qpart+1.0)/qparts * ddN.DT;
  }
  for(i=startdoc; i<enddoc; i++) {
    int  thisw =  add_doc(i, GibbsNone);
    int  r;
    if ( thisw<=1 ) {
      remove_doc(i, GibbsNone);
      continue;
    }
    if ( ddP.bdk!=NULL ) 
      misi_build(&dD, i, 0);
    map_query(i, mimap, found_buf);
    for (j=0; j<ddP.n_words; j++) {
      topcnt_buf[j] = 0;
      topwordscore_buf[j] = 0;
    }
    
    for (r=0; r<ddP.queryiter; r++) {
      gibbs_lda(GibbsNone, ddN.T, i, ddD.NdT[i], fact, &dD, 0, 0);
      query_docprob(i, mimap, fact, &dD, topcnt_buf, topwordscore_buf);
    }  
    /*
     *  now adjust stats
     */
    for (j=0; j<ddP.n_query; j++) 
      logprob[j] = 0;
    for (j=0; j<ddP.n_words; j++) {
      if ( wordunused[j]>0 )
	continue;
      if ( ddP.query[ddP.qword[j]]==j ) {
	topcnt_buf[j] /= ddP.queryiter;
	topwordscore_buf[j] /= ddP.queryiter;
      } else {
	/*  word in previous query so copy  */
	int jj =  ddP.query[ddP.qword[j]];
	topcnt_buf[j] = topcnt_buf[jj];
	topwordscore_buf[j] = topwordscore_buf[jj];
	found_buf[j] = found_buf[jj];
      }
      if ( wordunused[j]==0 )
	logprob[ddP.qid[j]] += topwordscore_buf[j];
    }
    if ( dots>0 && i>0 && (i%dots==0) ) 
      yap_message(".");
    if ( ddP.bdk!=NULL ) misi_unbuild(&dD,i,0);
    remove_doc(i, GibbsNone);
    /*
     *   enter into the arrays
     */
    for (j=0; j<ddP.n_query; j++) {
      if ( i<K || logprob[j] < topscore[j*K+topind[j*K+K-1]] ) {
	int newind, l;
	/*
	 *   better than current lowest 
	 */
	newind = bubble((i<K)?(i+1):K, 
			&topind[j*K], &topscore[j*K], logprob[j]);
	/*
	 *   save the current details
	 */
	topscore[j*K+newind] = logprob[j];
	topk[j*K+newind] = i;
	for (l=ddP.qposn[j]; l<ddP.qposn[j+1]; l++) {
	  topcnt[newind*ddP.n_words+l] = topcnt_buf[l]; 
	  topwordscore[newind*ddP.n_words+l] = topwordscore_buf[l]; 
	  found[newind*ddP.n_words+l] = found_buf[l]; 
	}
      }
    }
  }
  if ( dots>0 ) yap_message("\n");
  
  /*
   *  load df
   */
  df = calloc(ddN.W,sizeof(df[0]));
  if ( !df ) 
    yap_quit("Cannot allocate memory in gibbs_query()\n");
  n_df = data_df(stem,df);
  
  /*
   *  write result
   */
  {
    float *ws = fvec(ddP.n_words);
    FILE *fp = fopen(qname,"w");
    int q;
    if ( !fp )
      yap_sysquit("Cannot write query results to '%s'\n", qname);
    for (q=0; q<ddP.n_query; q++) {
      int nw = ddP.qposn[q+1]-ddP.qposn[q];
      for (i=0; i<K && i<ddN.DT && topk[topind[q*K+i]]>=0; i++) {
	int l, ind = topind[q*K+i];
	double tfidf;
        tfidf = bm25(topk[q*K+ind],&found[ind*ddP.n_words+ddP.qposn[q]],
			    &ddP.qword[ddP.qposn[q]], nw, ws);
	assert(ind>=0 && ind<K);
	fprintf(fp, "%d %d ", q, topk[q*K+ind]);
	fprintf(fp, "%.4f %.4lf ", topscore[q*K+ind]/nw, tfidf);
        if ( verbose>1 ) {
          for (l=ddP.qposn[q]; l<ddP.qposn[q+1]; l++)
            fprintf(fp, "%d ", found[ind*ddP.n_words+l]);
          for (l=ddP.qposn[q]; l<ddP.qposn[q+1]; l++)
            fprintf(fp, "%f ", topcnt[ind*ddP.n_words+l]);
          for (l=ddP.qposn[q]; l<ddP.qposn[q+1]; l++)
            fprintf(fp, "%f ", topwordscore[ind*ddP.n_words+l]);
          for (l=0; l<nw; l++)
            fprintf(fp, "%lf ", ws[l]);
        }
        fprintf(fp, "\n");
      }
    }
    fclose(fp);
    free(ws);
  }

  /*
   *  clean up
   */
  free(df);
  free(fact);
  if ( ddP.bdk!=NULL ) misi_free(&dD);
  if ( mimap ) free(mimap);
  free(found);
  free(topwordscore);
  free(topcnt);
  free(found_buf);
  free(topwordscore_buf);
  free(topcnt_buf);
  free(topscore);
  free(topind);
  free(topk);
  free(logprob);
}


