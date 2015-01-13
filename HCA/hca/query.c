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
#ifdef H_THREADS
#include <pthread.h>
#endif

#include "yap.h"
#include "util.h"
#include "fvec.h"
#include "hca.h"
#include "data.h"
#include "pctl.h"
#include "probs.h"
#include "termstats.h"
#ifdef H_THREADS
#include <pthread.h>
#endif
#include "atomic.h"

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
 *   return index where the new entry is placed;
 *   update topind[]
 *   K is length of vectors after insertion done;
 */
static int bubble(int K, int *topind, float *score, float newscore) {
  int k;
  int newind;
  if ( K>1 ) {
    int scale;
#if 0
    k = K-1;
    if ( !finite(score[topind[k]]) )
      k--;
    yap_message("bubble-in:");
    for (; k>0; k--) {
      yap_message("(%d)%f ", topind[k], score[topind[k]]);
      assert(score[topind[k-1]]<=score[topind[k]]);
    }
    yap_message("\n");
#endif
    /*  initial bisection search */
    newind = topind[K-1];
    assert(newscore<score[newind]);
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
  } else {
    newind = 0;
    k = 0;
  }
  topind[k] = newind;
#if 0
  score[newind] = newscore;
    yap_message("bubble-out:");
  for (k=K-1; k>0; k--) {
    yap_message("(%d)%f ", topind[k], score[topind[k]]);
    assert(score[topind[k-1]]<=score[topind[k]]);
  }
  yap_message("\n");
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

/*
 *  build query map for document d
 */
static void map_query(int d, int *map, int *found) {
  int mi = 0;
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

/************************************************************************
 *    temporary versions of results for when gibbs running on single doc
 */
typedef struct QD_s {
  int     *found;
  float   *topcnt;
  float   *topwordscore;
  double  *logprob;
  /*
   *    mapping from query word posn. to its mi in current doc
   *       >ddN.N  = not in current doc
   *       -ve  = has no mi since occurs just once, found at
   *              posn  (-map[]-1)
   *       non -ve = mi value
   */
  int     *map;
} QD_t;
static void QD_init(QD_t *buf) {
  if ( ddP.bdk!=NULL ) 
    buf->map = malloc(sizeof(buf->map[0])*ddP.n_words);
  else
    buf->map = NULL;
  buf->topcnt = malloc(sizeof(buf->topcnt[0])*ddP.n_words);
  buf->topwordscore = malloc(sizeof(buf->topwordscore[0])*ddP.n_words);
  buf->found = malloc(sizeof(buf->found)*ddP.n_words);
  buf->logprob = malloc(sizeof(buf->logprob[0])*ddP.n_query);
  if ( !buf->topcnt || !buf->topwordscore || !buf->found || !buf->logprob )
    yap_quit("Cannot allocate memory in gibbs_query()\n");
}
static void QD_free(QD_t *buf) {
  free(buf->found);
  free(buf->topwordscore);
  free(buf->topcnt);
  free(buf->logprob);
  if ( buf->map ) free(buf->map);
}

/*
 *  single document query processing
 */
static void query_run(int i, QD_t *buf, D_MiSi_t *dD, char *wordunused) {
  int  thisw =  add_doc(i, GibbsNone);
  int  r, j;
  float *fact = fvec(ddN.T*4);

  if ( thisw<=1 ) {
    remove_doc(i, GibbsNone);
    return;
  }
  if ( ddP.bdk!=NULL ) 
    misi_build(dD, i, 0);
  map_query(i, buf->map, buf->found);
  for (j=0; j<ddP.n_words; j++) {
    buf->topcnt[j] = 0;
    buf->topwordscore[j] = 0;
  }
    
  for (r=0; r<ddP.queryiter; r++) {
    gibbs_lda(GibbsNone, ddN.T, i, ddD.NdT[i], fact, dD, 0, 0);
    query_docprob(i, buf->map, fact, dD, buf->topcnt, buf->topwordscore);
  }  
  /*
   *  now adjust stats
   */
  for (j=0; j<ddP.n_query; j++) 
    buf->logprob[j] = 0;
  for (j=0; j<ddP.n_words; j++) {
    if ( wordunused[j]>0 )
      continue;
    if ( ddP.query[ddP.qword[j]]==j ) {
      buf->topcnt[j] /= ddP.queryiter;
      buf->topwordscore[j] /= ddP.queryiter;
    } else {
      /*  word in previous query so copy  */
      int jj =  ddP.query[ddP.qword[j]];
      buf->topcnt[j] = buf->topcnt[jj];
      buf->topwordscore[j] = buf->topwordscore[jj];
      buf->found[j] = buf->found[jj];
    }
    if ( wordunused[j]==0 )
      buf->logprob[ddP.qid[j]] += buf->topwordscore[j];
  }
  if ( ddP.bdk!=NULL ) misi_unbuild(dD,i,0);
  remove_doc(i, GibbsNone);  
  free(fact);
}

/********************************************************
 *    top data so far for query
 *    these store statistics of the results, for printing
 *    these are unordered, ordered by ind[]
 */
typedef struct QT_s {
  /*
   *  document details for top results
   */
  /*      results saved (per query) */
  int     *saved;
  /*      document score  */
  float   *score;
  /*      document number   */
  int     *k;
  /*      an index into k[] which maintains ordering  */
  int     *ind;  
  /*
   *    per word stats for top results saved
   */
  int     *found;
  float   *cnt;
  float   *wordscore;
  /*      flags if word is irrelevant, thus not scored  */
  char    *wordunused;
  /*
   *      locks to add values;
   */
#ifdef H_THREADS
  pthread_mutex_t *mutex;
#endif

} QT_t;
static void QT_init(QT_t *buf, int topQ, int procs) {
  int i;
  buf->wordunused = malloc(sizeof(buf->wordunused[0])*ddP.n_words);
  buf->cnt = malloc(sizeof(buf->cnt[0])*topQ*ddP.n_words);
  buf->wordscore = malloc(sizeof(buf->wordscore[0])*topQ*ddP.n_words);
  buf->score = malloc(sizeof(buf->score[0])*topQ*ddP.n_query);
  buf->saved = malloc(sizeof(buf->saved[0])*ddP.n_query);
  buf->ind = malloc(sizeof(buf->ind[0])*topQ*ddP.n_query);
  buf->k = malloc(sizeof(buf->k[0])*topQ*ddP.n_query);
  buf->found = malloc(sizeof(buf->found)*ddP.n_words*topQ);
  if ( !buf->cnt || !buf->wordscore || !buf->found || !buf->ind 
       ||  !buf->k ||  !buf->score ||  !buf->wordunused || !buf->saved )
    yap_quit("Cannot allocate memory in gibbs_query()\n");
  for (i=0; i<ddP.n_words; i++) {
    buf->wordunused[i] = 0;
  }
  for (i=0; i<ddP.n_query; i++) {
    buf->saved[i] = 1;
  }
  for (i=0; i<topQ*ddP.n_query; i++) {
    buf->ind[i] = i%topQ;
    buf->k[i] = -1;
    buf->score[i] = HUGE_VAL;
  }
#ifdef H_THREADS
  if ( procs>1 ) {
    buf->mutex = malloc(sizeof(buf->mutex[0])*ddP.n_query);
    if ( !buf->mutex )
      yap_quit("Cannot allocate memory in gibbs_query()\n");
    for (i=0; i<ddP.n_query; i++) {
      pthread_mutex_init(&buf->mutex[i], NULL);
    }
  } else
    buf->mutex = NULL;
#endif
}
static void QT_free(QT_t *buf) {
  free(buf->wordscore);
  free(buf->cnt);
  free(buf->saved);
  free(buf->score);
  free(buf->ind);
  free(buf->k);
  free(buf->found);
  free(buf->wordunused);
#ifdef H_THREADS
  if ( buf->mutex ) {
    int i;
    for (i=0; i<ddP.n_query; i++) {
      pthread_mutex_destroy(&buf->mutex[i]);
    }
    free(buf->mutex);
  }
#endif
}
static void QT_save(int i, int topQ, QT_t *top, QD_t *doc) {
  int j;
  /*
   *   enter into the arrays
   */
  for (j=0; j<ddP.n_query; j++) {
    /*
     *   top->saved[j] is non-decreasing and
     *   top->score[j*topQ+top->ind[j*topQ+topQ-1]] is non-increasing
     *   so we test once outside of lock, and retest inside lock
     */   
    if ( top->saved[j]<topQ || 
	 doc->logprob[j] < top->score[j*topQ+top->ind[j*topQ+topQ-1]] ) {
      int newind, l;
#ifdef H_THREADS
      if ( top->mutex ) {
	pthread_mutex_lock(&top->mutex[j]);
	if ( top->saved[j]>=topQ &&
	     doc->logprob[j]>=top->score[j*topQ+top->ind[j*topQ+topQ-1]] ) {
	  /*
	   *  so above test result changed in new lock state, so quit early
	   */
	  pthread_mutex_unlock(&top->mutex[j]);
	  continue;
	}
      }
#endif
      /*
       *   must insert
       */
      newind = bubble(top->saved[j],
		      &top->ind[j*topQ], &top->score[j*topQ], doc->logprob[j]);
      /*
       *   save the current details
       *     (all specific to query j)
       */
      top->score[j*topQ+newind] = doc->logprob[j];
      top->k[j*topQ+newind] = i;
      for (l=ddP.qposn[j]; l<ddP.qposn[j+1]; l++) {
	// 	yap_message("INSERT: n=%d l=%d index=%d\n", newind, l, newind*ddP.n_words+l);
	  top->cnt[newind*ddP.n_words+l] = doc->topcnt[l]; 
	  top->wordscore[newind*ddP.n_words+l] = doc->topwordscore[l]; 
	  top->found[newind*ddP.n_words+l] = doc->found[l]; 
      }
      if ( top->saved[j]<topQ )
	top->saved[j]++;
#ifdef H_THREADS
      if ( top->mutex ) {
	pthread_mutex_unlock(&top->mutex[j]);
      }
#endif
    }
  }
}

static void pprune(double *pvec, int N, float cutoff) {
  int i;
  double maxp = 0;
  for (i=0; i<N; i++)
    if ( pvec[i]>maxp )
      maxp = pvec[i];
  for (i=0; i<N; i++)
    if ( pvec[i]<maxp*cutoff )
      pvec[i] = 0;
}

/*
 *   compute doc probability vector for query q
 */
static float *qdocprob(int topQ, QT_t *top, int q) {
  float *dp = fvec(topQ);
  int i;
  double minscore, totscore = 0;
  int *docs = &top->k[q*topQ];
  float *score = &top->score[q*topQ];

  minscore = score[0];
  for (i=0; i<topQ; i++) 
    if ( docs[i]>=0 && minscore>score[i] )
      minscore = score[i];
  for (i=0; i<topQ; i++) {
    if ( docs[i]<0 )
      continue;
    totscore += dp[i] = exp(-(score[i]-minscore));
    //yap_message(" %f %f", score[i], dp[i]);
  }
  if ( totscore==0 ) {
    free(dp);
    return NULL;
  }
  //yap_message("\n");
  for (i=0; i<topQ; i++) {
    if ( docs[i]>=0 ) {
      dp[i] /= totscore;
      //yap_message(" %f", dp[i]);
    }
  }
  //yap_message("\n");
  return dp;
}

/*
 *  write result
 */
static void QT_write(char *qname, int topQ, QT_t *top) {
  float *ws = fvec(ddP.n_words);
  FILE *fp = fopen(qname,"w");
  int i, q;
  if ( !fp )
    yap_sysquit("Cannot write query results to '%s'\n", qname);
  for (q=0; q<ddP.n_query; q++) {
    int nw = ddP.qposn[q+1]-ddP.qposn[q];
    float *dp = qdocprob(topQ, top, q);
    if ( !dp )
      yap_quit("No results for query %d\n", q);
    yap_message("Effective #docs for query %d = %lf\n", q+1, 
		exp(fv_entropy(dp,topQ)));
    for (i=0; i<top->saved[q]; i++) {
      int l, ind = top->ind[q*topQ+i];
      double tfidf;
      if ( top->k[q*topQ+ind]<0 )
	continue;
      tfidf = bm25(top->k[q*topQ+ind],
		   &top->found[ind*ddP.n_words+ddP.qposn[q]],
		   &ddP.qword[ddP.qposn[q]], nw, ws);
      assert(ind>=0 && ind<topQ);
      fprintf(fp, "%d %d ", q, top->k[q*topQ+ind]);
      fprintf(fp, "%f %.4f %.4lf ", dp[ind], top->score[q*topQ+ind]/nw, tfidf);
      if ( verbose>1 ) {
	for (l=ddP.qposn[q]; l<ddP.qposn[q+1]; l++)
	  fprintf(fp, "%d ", top->found[ind*ddP.n_words+l]);
	for (l=ddP.qposn[q]; l<ddP.qposn[q+1]; l++)
	  fprintf(fp, "%f ", top->cnt[ind*ddP.n_words+l]);
	for (l=ddP.qposn[q]; l<ddP.qposn[q+1]; l++)
	  fprintf(fp, "%f ", top->wordscore[ind*ddP.n_words+l]);
	for (l=0; l<nw; l++)
	  fprintf(fp, "%lf ", ws[l]);
      }
      fprintf(fp, "\n");
    }
    free(dp);
  }
  fclose(fp);
  free(ws);
}

/*
 *  write term report;
 *  rather inefficient use of memory
 */
#define TERM_CUTOFF 1.0e-4
static void QT_terms(char *stem, char *qname, int topQ, QT_t *top) {
  FILE *fp = fopen(qname,"w");
  int i, q;
  double *dprob;
  double *wprob, *tprob;
  if ( !fp )
    yap_sysquit("Cannot write query term data to '%s'\n", qname);

  wprob = malloc(sizeof(wprob[0])*ddN.W);
  dprob = malloc(sizeof(dprob[0])*ddN.DT);
  tprob = malloc(sizeof(tprob[0])*ddN.T);
  if ( !dprob || !wprob || !tprob )
    yap_quit("Out of memory in QT_terms()\n");

  for (q=0; q<ddP.n_query; q++) {
    int *docs = &top->k[q*topQ];
    float *dp = qdocprob(topQ, top, q);
    T_stats_t *tstats;
    /*
     *  build set of doc probabilities in dprob
     */
    
    for (i=0; i<ddN.DT; i++) 
      dprob[i] = 0;
    for (i=0; i<topQ; i++) {
      if ( docs[i]<0 )
	continue;
      dprob[docs[i]] = dp[i];
    }
    /*
     *  build word probs
     */
    for (i=0; i<ddN.W; i++) 
      wprob[i] = 0;
    for (i=0; i<ddN.DT; i++) {
      int l;
      if ( dprob[i]==0 )
	continue;
      for (l=ddD.NdTcum[i]; l<ddD.NdTcum[i+1]; l++) {
	wprob[ddD.w[l]] += dprob[i];
      }
    }
    pprune(wprob,ddN.W,TERM_CUTOFF);
    for (i=0; i<ddN.W; i++) 
      if ( wprob[i]>0 ) {
	double wp = wprob[i]/top->saved[q];
	fprintf(fp, "W %d %d %lf", q, i, wp);
	if ( ddN.tokens )
	  fprintf(fp, " %s", ddN.tokens[i]);
	fprintf(fp, "\n");
      }
    /*
     *  build term probs
     */
    tstats = twstats_init(dprob, ddD.NdTcum, ddN.T, ddN.DT, stem);
    pprune(tstats->Nt,tstats->K,TERM_CUTOFF);
    for (i=0; i<tstats->K; i++) 
      if ( tstats->Nt[i]>0 ) {
	double wp = tstats->Nt[i]/top->saved[q];
	fprintf(fp, "K %d %d %lf", q, i, wp);
	fprintf(fp, " %s", tstats->tokens[i]);
	fprintf(fp, "\n");
      }
    tstats_free(tstats);
    /*
     *  build topic probs
     */
    for (i=0; i<ddN.T; i++) 
      tprob[i] = 0;
    for (i=0; i<ddN.DT; i++) {
      int Nt = ddD.NdT[i];
      int l;
      if ( dprob[i]==0 )
	continue;
      for (l=ddD.NdTcum[i]; l<ddD.NdTcum[i+1]; l++) {
	tprob[Z_t(ddS.z[l])] += dprob[i]/Nt;
      }
    }
    pprune(tprob,ddN.T,TERM_CUTOFF);
    for (i=0; i<ddN.T; i++) 
      if ( tprob[i]>0 ) {
	fprintf(fp, "T %d %d %lf\n", q, i, tprob[i]);
      }
  }
  fclose(fp);
  free(wprob);
  free(dprob);
  free(tprob);
}

/*
 *   arguments for parallel call to Gibbs querying sampler
 */
typedef struct D_qargs_s {
  int dots;
  int processid;
  int procs;
  int *doc;
  QT_t *top;  /*  top results so far */
  int topQ;
} D_qargs_p;

void *querying_p(void *qargs) {
  int i;
  D_MiSi_t dD;  
  D_qargs_p *par =(D_qargs_p *) qargs;
  /*  results from single query */
  QD_t QDbuf;

  if ( ddP.bdk!=NULL ) misi_init(&ddM,&dD);
  QD_init(&QDbuf);

  while ( (i=atomic_add(*par->doc,1))<ddN.DT ) {   
    query_run(i, &QDbuf, &dD, par->top->wordunused);
    if ( par->dots>0 && i>0 && (i%par->dots==0) ) 
      yap_message(".");
    QT_save(i, par->topQ, par->top, &QDbuf);
  }

  QD_free(&QDbuf);
  if ( ddP.bdk!=NULL ) misi_free(&dD);
  return NULL;
}

/*
 *    run regular gibbs cycles on the data with phi used;
 *    the evaluation on each doc, and sample word probs
 *
 *    topQ = number of top results to retain
 */
void gibbs_query(char *stem, int topQ, char *queryfile, int dots, int procs,
		 int excludedoc, float scale) {
  /*  top results so far */
  QT_t Qtop;
#ifdef H_THREADS
  pthread_t thread[procs];
#endif
  D_qargs_p parg[procs];
  int doc;
  int pro;
  char *qdname = yap_makename(queryfile,".docs");
  char *qtname = yap_makename(queryfile,".terms");

  QT_init(&Qtop, topQ, procs);
  
  /*
   *  check words to exclude using topics
   */
  if ( ddP.n_excludetopic>0 ) {
    int i;
    double *tprob = malloc(sizeof(tprob[0])*ddN.T);
    assert(ddS.Ndt);
    get_probs(tprob);
    yap_probs();
    if ( verbose>1 )
      yap_message("Excluding words: ");
    for (i=0; i<ddP.n_words; i++) {
      int t = besttopic(ddP.qword[i],tprob);
      if ( Q_excludetopic(t) ) {
	Qtop.wordunused[i] = 1;
	if ( verbose>1 )
	  yap_message(" %d/%d", (int)ddP.qword[i], t);
      }
    } 
    if ( verbose>1 )
      yap_message("\n");
    free(tprob);
  }
  
  doc = -1;
  for (pro = 0 ; pro < procs ; pro++){
    parg[pro].dots=dots;
    parg[pro].processid=pro;
    parg[pro].procs=procs;
    parg[pro].doc = &doc;
    parg[pro].top = &Qtop;
    parg[pro].topQ = topQ;
#ifndef H_THREADS
    querying_p(&parg[pro]);
#else
    if ( procs==1 )
      querying_p(&parg[pro]);
    else if ( pthread_create(&thread[pro],NULL,querying_p,(void*) &parg[pro]) != 0){
      yap_message("thread failed %d\n",pro+1 );
    }
#endif
  }
#ifdef H_THREADS
  if ( procs>1 ) {
    //waiting for threads to finish
    for (pro = 0; pro < procs; pro++){
      pthread_join(thread[pro], NULL);
    }
  }
#endif

  if ( excludedoc ) {
    int q, i;
    for (q=0; q<ddP.n_query; q++) {
      for (i=0; i<Qtop.saved[q] && i<excludedoc; i++) {
	int ind = Qtop.ind[q*topQ+i];
	//  delete i-th entry
	Qtop.k[q*topQ+ind] = -1;
	Qtop.score[q*topQ+ind] = HUGE_VAL;
      }
      //  shuffle down excludedoc entries
      for (i=0; i<Qtop.saved[q]-excludedoc; i++) {
	Qtop.ind[q*topQ+i] = Qtop.ind[q*topQ+i+excludedoc];
      }
    }
  }
  if ( scale != 1.0 ) {
    int i;
    for (i=0; i<topQ*ddP.n_query; i++)
      Qtop.score[i] *= scale;
  }

  if ( dots>0 ) yap_message("\n");
  
  /*
   *  load df
   */
  df = calloc(ddN.W,sizeof(df[0]));
  if ( !df ) 
    yap_quit("Cannot allocate memory in gibbs_query()\n");
  n_df = data_df(stem, df);
  
  QT_write(qdname, topQ, &Qtop);
  QT_terms(stem, qtname, topQ, &Qtop);

  /*
   *  clean up
   */
  free(df);
  free(qdname);
  free(qtname);
  QT_free(&Qtop);
}


