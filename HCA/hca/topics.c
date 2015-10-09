/*
 * Auxiliary topic routines for display
 * Copyright (C) 2012-2013 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *
 *   First tried scoring terms based on count.
 *   Then modified things to an expected value score,
 *   but applied a filter to remove low predictive terms.
 *     
 */
/*  KL  */
#define KL

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <assert.h>
#include "hca.h"
#include "stats.h"
#include "lgamma.h" 
#include "util.h" 
#include "yap.h" 
#include "diag.h" 
#include "probs.h" 
#include "pmi.h" 
#include "fvec.h" 
#include "termstats.h" 

/*************************************************
 *  return sort order for topics by size
 */
static float *stvec;
static int pcompar(const void *a, const void *b) {
  float na = stvec[*(uint32_t*)a];
  float nb = stvec[*(uint32_t*)b];
  if ( na<nb ) return 1;
  if ( na>nb ) return -1;
  return 0;
}
static uint32_t *sorttops(float *vec, int K) {
  uint32_t *psort = u32vec(K);
  int k;
  for (k=0; k<K; k++) psort[k] = k;
  stvec = vec;
  qsort(psort, K, sizeof(*psort), pcompar);
  return psort;
}

/*************************************************
 *   build topic proportions vec for doc d
 */
static float *topprop(int d) {
  float *vec = fvec(ddN.T);
  double tot = 0;
  int k;
  if ( !vec) return NULL;
  if ( ddP.theta ) {
    for (k=0; k<ddN.T; k++) {
      assert( ddP.theta[d][k]>=0 );
      tot += vec[k] = ddP.theta[d][k];
    }
  } else { 
    for (k=0; k<ddN.T; k++) {
      assert(ddS.Ndt[d][k]>=0);
      tot += vec[k] = ddS.Ndt[d][k];
    }
  }
  if ( tot <=0 )
    return vec;
  for (k=0; k<ddN.T; k++)
    vec[k] /= tot;
  return vec;
}

/*
 *   global probability for topic:
 *   taken from data if exists;
 *   else from ddP.alphapr
 *   else is 0
 */
static float *globalprop() {
  float *vec = fvec(ddN.T);
  double tot = 0;
  int k;
  if ( !vec) return NULL;
  if ( ddS.Nwt  ) {
    for (k=0; k<ddN.T; k++) {
      tot += vec[k] = ddS.NWt[k];
    }
  } else if ( ddS.Ndt ) {
    /*  have to total from Ndt */
    int d;
    uint32_t *uvec = u32vec(ddN.T);
    for (d=0; d<ddN.DT; d++) {
      for (k=0; k<ddN.T; k++) 
	uvec[k] += ddS.Ndt[d][k];
    }
    for (k=0; k<ddN.T; k++) 
      tot += vec[k] = uvec[k];
    free(uvec);
  } else if ( ddP.alphapr ) {
    for (k=0; k<ddN.T; k++) 
      tot += vec[k] = ddP.alphapr[k];
  }
  if ( tot <=0 )
    return vec;
  for (k=0; k<ddN.T; k++) 
    vec[k] /= tot;
  return vec;
}

/*
 *   Build a matrix giving correlation of topic proportions for docs
 *   This is O(D*T*T) so very slow.
 */
float **hca_topmtx() {
  float **mtx = fmat(ddN.T, ddN.T);
  float *vec = fvec(ddN.T);
  int i, t1, t2;
  for (i=0; i<ddN.DT; i++) {
    float *tvec = topprop(i);
    for (t1=0; t1<ddN.T; t1++) 
      vec[t1] += tvec[t1];
    for (t1=0; t1<ddN.T; t1++) {
      if ( tvec[t1]<=0 )
        continue;
      for (t2=0; t2<t1; t2++) 
	mtx[t1][t2] += tvec[t1] * tvec[t2];
    }
    free(tvec);
  }
  for (t1=0; t1<ddN.T; t1++)
    vec[t1] /= ddN.DT;
  for (t1=0; t1<ddN.T; t1++)
    for (t2=0; t2<t1; t2++) {
      mtx[t1][t2] = mtx[t1][t2]/(ddN.DT*vec[t1]*vec[t2]);
    }
  free(vec);
  return mtx;
}

/*
 *   Build a vector giving number of times each topic is the
 *   most common in a doc.
 */
static uint32_t *hca_top1cnt() {
  uint32_t *cnt = u32vec(ddN.T);
  int i, t;
  for (i=0; i<ddN.DT; i++) {
    float *tvec = topprop(i);
    int maxt = 0;
    for (t=0; t<ddN.T; t++) 
      if ( tvec[t]>tvec[maxt] )
	maxt = t;
    cnt[maxt]++;
    free(tvec);
  }     
  return cnt;
}

/*************************************************
 *   build document proportions vec for topic k
 */
static float *docprop(int k) {
  float *vec = fvec(ddN.DT);
  int i;
  if ( !vec) return NULL;
  if ( ddP.theta ) {
    for (i=0; i<ddN.DT; i++)
      vec[i] = ddP.theta[i][k];
  } else { 
    for (i=0; i<ddN.DT; i++) {
      if ( ddS.NdT[i]>0 )
        vec[i] = ddS.Ndt[i][k] / (float)ddS.NdT[i];
      else 
        vec[i] = 0;
    }
  }
  return vec;
}
/*************************************************
 *  find top-K by bubble sort
 */
static void topk(int K, int N, uint32_t *ind, double (*foo)(int)) {
  int k, n;
  double val;
  for (k=0; k<K && k<N-1; k++) {
    val = foo(ind[k]);
    for (n=k+1; n<N; n++)
      if ( foo(ind[n])>val ) {
	int ss = ind[k];
	ind[k] = ind[n];
	ind[n] = ss;
	val = foo(ind[k]);
      }
  }
}


/**************************************************
 *  statistics of word counts in data
 */
static uint32_t NWK = 0;
static uint32_t *NwK = NULL;
static uint32_t termNWK = 0;
static uint32_t *termNwK = NULL;
static uint32_t **termNwk = NULL;

static void build_NwK() {
  int w, k;
  NwK = u32vec(ddN.W);
  if ( !NwK )
    yap_quit("Out of memory in hca_displaytopics()\n");
  for (w=0; w<ddN.W; w++) {
    NwK[w] = 0;
  }
  NWK = 0;
  if ( !ddS.Nwt ) {
    /*
     *  recompute from scratch
     */
    int i;
    for (i=0; i<ddN.NT; i++) 
      NwK[ddD.w[i]]++;
    for (w=0; w<ddN.W; w++) 
      NWK += NwK[w];
  } else {
    for (w=0; w<ddN.W; w++) {
      for (k=0; k<ddN.T; k++) {
	NwK[w] += ddS.Nwt[w][k];    //  should use CCT_ReadN()
      }
      NWK += NwK[w];
    }
  }
  if ( NWK==0 )
    yap_quit("empty NWK in build_NwK()\n");
}
static void build_termNwK(T_stats_t *ptr) {
  int w, k;
  termNwK = u32vec(ptr->K);
  if ( !termNwK )
    yap_quit("Out of memory in hca_displaytopics()\n");
  for (w=0; w<ptr->K; w++) {
    termNwK[w] = 0;
  }
  termNWK = 0;
  for (w=0; w<ptr->K; w++) {
    for (k=0; k<ddN.T; k++) {
      termNwK[w] += ptr->Nkt[w][k]; 
    }
    termNWK += termNwK[w];
  }
  if ( termNWK==0 ) 
    yap_quit("empty termNWK in build_termNwK(), collocations empty!\n");
  termNwk = ptr->Nkt;
}

/***************************************************
 *
 *   scoring utilities
 */

/*
 *  all the scoring functions below require this static set to work
 */
static int tscorek;

static unsigned getn(int w) {
  if ( tscorek<0 ) 
    return ddS.TwT[w];
  else
    return ddS.Nwt[w][tscorek];
}
static unsigned termgetn(int w) {
  return termNwk[w][tscorek];
}

static double idfscore(int w) {
  if ( tscorek>=0 && ddP.phi ) 
    return ddP.phi[tscorek][w]/NwK[w];
  return (getn(w)+0.2)/(NwK[w]+0.2*ddN.T);
}
static double termidfscore(int w) {
  return (termgetn(w)+0.2)/(termNwK[w]+0.2*ddN.T);
}

static double phiscore(int w) {
  if ( tscorek<0 ) {
    if ( ddP.betapr ) 
      return ddP.betapr[w];
    return ddS.TwT[w];
  }
  assert(ddP.phi || ddS.phi);
  if ( ddP.phi )
    return ddP.phi[tscorek][w];
  return ddS.phi[tscorek][w];
}

static double phiratioscore(int w) {
  assert(ddS.TwT);
  assert(tscorek>=0);
  return wordprob(w,tscorek)/ddS.TwT[w];
}

static double lowerQ;
static double Qscore(int w) {
  int n = getn(w);
  double N;
  N = NwK[w];
  if ( n/((double)N) <= lowerQ )
    return 0;
  return N/((double)ddN.NT) * 
    (n/((double)N)-lowerQ)/(1-lowerQ);
}

static double costscore(int w) {
  return getn(w)/((double)ddS.NWt[tscorek]) 
    - (NwK[w]*ddS.NWt[tscorek]/((double)ddN.NT))/((double)ddN.NT);
}

static double countscore(int w) {
  return getn(w);
}
static double termcountscore(int w) {
  return termgetn(w);
}

/******************************************************
 *  coherence calcs
 */
double coherence(uint32_t **mtx, int N) {
  int i1, i2; 
  double co = 0;
  for (i1=0; i1<N; i1++)
    for (i2=0; i2<i1; i2++)
      co += log((mtx[i1][i2]+1.0)/mtx[i1][i1]);
  co /= N*(N-1)/2;
  return co;
}

double coherence_word(uint32_t **mtx, int N, int i1) {
  int i2; 
  double co = 0;
  if (i1==0 )
    return 0.0;
  assert(i1<N);
  for (i2=0; i2<i1; i2++)
    co += log((mtx[i1][i2]+1.0)/mtx[i1][i1]);
  co /= i1;
  return co;
}

/********************************************************
 *  build non-zero indices for topic
 */
static int buildindk(int k, uint32_t *indk) {
  int w;
  int cnt=0;
  if ( k<0 ) {
    /*
     *  for root topic
     */
    if ( ddP.betapr ) {
      for (w=0; w<ddN.W; w++) 
	if ( ddP.betapr[w]>0.05/ddN.W ) indk[cnt++] = w;
      return cnt;
    }
    if ( ddP.phi!=NULL ) 
      return 0;
    assert(ddS.TwT);
    for (w=0; w<ddN.W; w++) {
      if ( ddS.TwT[w]>0 ) indk[cnt++] = w;
    }
    return cnt;
  }
  if ( ddP.phi ) {
    for (w=0; w<ddN.W; w++) 
      if ( ddP.phi[k][w]>0.05/ddN.W ) indk[cnt++] = w;
  } else {
    assert(ddS.Nwt);
    for (w=0; w<ddN.W; w++) {
      if ( ddS.Nwt[w][k]>0 ) indk[cnt++] = w;
    }
  } 
  return cnt;
}
static int buildtermindk(int t, uint32_t *indk, T_stats_t *ts) {
  int w;
  int cnt=0;
  for (w=0; w<ts->K; w++) 
    if ( ts->Nkt[w][t]>0 ) indk[cnt++] = w;
  return cnt;
}

uint32_t **classbytopic(char *resstem) {
  int i;
  uint32_t **TbyC;
  /*
   *  write topic by class confusion matrix
   */
  TbyC = u32mat(ddN.T,ddN.C);
  for (i=0; i<ddN.NT; i++) {
    int d = ddD.d[i];
    TbyC[Z_t(ddS.z[i])][ddD.c[d]]++;
  }
  if ( resstem ) {
    char *fname;
    fname = yap_makename(resstem,".tbyc");
    write_u32sparse(ddN.T,ddN.C,TbyC,fname);
    free(fname);
  }
  return TbyC;
}

void hca_displayclass(char *resstem) {
  int i, k;
  double ent = 0;
  uint32_t **TbyC;

  /*
   *  write topic by class confusion matrix
   */
  TbyC = classbytopic(resstem);

  /*
   *  now report entropies
   */
  yap_message("Class entropies by topic: ");
  for (k=0; k<ddN.T; k++) {
    double me = 0;
    double tot = 0;
    for (i=0; i<ddN.C; i++) 
      tot += TbyC[k][i];
    for (i=0; i<ddN.C; i++) {
      double p;
      if ( TbyC[k][i]>0 ) {
	p = ((double)TbyC[k][i])/tot;
	me -= p * log(p) * M_LOG2E;
      }
    }
    ent += me;
    yap_message(" %.3lf", me);
  }
  yap_message(" -> %.3lf\n", ent/ddN.T);
  free(TbyC[0]); free(TbyC);
}

void hca_displaytopics(char *stem, char *resstem, int topword, 
                       enum ScoreType scoretype, int pmicount, int fullreport) {
  int w,k;
  uint32_t *termindk = NULL;
  uint32_t *indk = NULL;
  int Nk_tot = 0;
  double (*termtscore)(int) = NULL;
  double (*tscore)(int) = NULL;
  double sparsityword = 0;
  double sparsitydoc = 0;
  double underused = 0;
  uint32_t *top1cnt = NULL;
  FILE *fp;
  float *tpmi = NULL;
  char *topfile;
  char *repfile;
  uint32_t *psort;
  FILE *rp = NULL;
  float *gtvec = globalprop();
  float *gpvec = calloc(ddN.W,sizeof(gpvec[0]));
  float *pvec = calloc(ddN.W,sizeof(pvec[0]));
#ifdef KL
  float *p0vec = calloc(ddN.W,sizeof(p0vec[0]));
#endif
  double *ngalpha = NULL;
  T_stats_t *termstats;
  
  if ( pmicount>topword )
    pmicount = topword;
  if ( scoretype == ST_idf ) {
    tscore = idfscore;
  } else if ( scoretype == ST_phirat ) {
    tscore = phiratioscore;
  } else if ( scoretype == ST_phi ) {
    tscore = phiscore;
  } else if ( scoretype == ST_count ) {
    tscore = countscore;
  } else if ( scoretype == ST_cost ) {
    tscore = costscore;
  } else if ( scoretype == ST_Q ) {
    tscore = Qscore;
    lowerQ = 1.0/ddN.T;
  }  

  if ( ddP.NGalpha ) {
    ngalpha = dvec(ddN.T);
    get_probs(ngalpha);
  }

  /*
   *  returns null if no relevant data file
   */
  termstats = tstats_init(ddS.z, ddD.NdTcum, ddN.T, ddN.DT, stem);
  if ( termstats ) {
    if ( scoretype == ST_idf ) {
      termtscore = termidfscore;
    } else 
      termtscore = termcountscore;
  }  

  
  /*
   *  first collect counts of each word/term,
   *  and build gpvec (mean word probs)
   */
  build_NwK();
  if ( termstats )
    build_termNwK(termstats);
  {
    /*
     *  gpvec[] is normalised NwK[]
     */
    double tot = 0;
    for (w=0; w<ddN.W; w++)
      tot += gpvec[w] = NwK[w]+0.1; 
    for (w=0; w<ddN.W; w++)
      gpvec[w] /= tot;
  }
  if ( ddS.Nwt ) {
    for (k=0; k<ddN.T; k++) {
      Nk_tot += ddS.NWt[k];
    }
  } 
  
  psort = sorttops(gtvec, ddN.T);
  
  top1cnt = hca_top1cnt();
  if ( !top1cnt )
    yap_quit("Cannot allocate top1cnt in hca_displaytopics()\n");

  if ( pmicount ) {
    tpmi = malloc(sizeof(*tpmi)*(ddN.T+1));
    if ( !tpmi )
      yap_quit("Cannot allocate tpmi in hca_displaytopics()\n");
  }
  indk = malloc(sizeof(*indk)*ddN.W);
  if ( !indk )
    yap_quit("Cannot allocate indk in hca_displaytopics()\n");
  if ( termstats ) {
    termindk = malloc(sizeof(*indk)*termstats->K);
    if ( !termindk )
      yap_quit("Cannot allocate termindk in hca_displaytopics()\n");
  }

#ifdef KL
  ddD.df = malloc(sizeof(*ddD.df)*ddN.W);
  ddD.n_df = data_df(stem,ddD.df);
  for (w=0; w<ddN.W; w++)
    p0vec[w] = ddD.df[w];
#endif
  
  /*
   *   two passes through, 
   *           first to build the top words and dump to file
   */
  repfile = yap_makename(resstem,".topset");
  topfile = yap_makename(resstem,".toplst");
  fp = fopen(topfile,"w");
  if ( !fp ) 
    yap_sysquit("Cannot open file '%s' for write\n", topfile);
  yap_message("\n");
  for (k=0; k<ddN.T; k++) {
    int cnt, termcnt = 0;
    tscorek = k;
    /*
     *    build sorted word list
     */
    cnt = buildindk(k, indk);
    topk(topword, cnt, indk, tscore);
    if ( cnt==0 )
      continue;
    if ( termstats ) {
      termcnt = buildtermindk(k, termindk, termstats);
      topk(topword, termcnt, termindk, termtscore);
    }
    /*
     *   dump words to file
     */
    fprintf(fp,"%d: ", k);
    for (w=0; w<topword && w<cnt; w++) {
      fprintf(fp," %d", (int)indk[w]);
    }
    if ( termstats ) {
      for (w=0; w<topword && w<termcnt; w++) {
	fprintf(fp," %d", (int)termstats->Kmin+termindk[w]);
      }
    }
    fprintf(fp, "\n");
  }
  if ( ddP.PYbeta && (ddP.phi==NULL || ddP.betapr)  ) {
    int cnt;
     /*
     *    dump root words
     */
    tscorek = -1;
    cnt = buildindk(-1, indk);
    topk(topword, cnt, indk, (ddP.phi==NULL)?countscore:phiscore);
    fprintf(fp,"-1:");
    for (w=0; w<topword && w<cnt; w++) {
      fprintf(fp," %d", (int)indk[w]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  if ( verbose>1 ) yap_message("\n");

  if ( pmicount ) {
    /*
     * compute PMI
     */
    char *toppmifile;
    char *pmifile;
    double *tp;
    tp = dvec(ddN.T);
    pmifile=yap_makename(stem,".pmi");
    toppmifile=yap_makename(resstem,".toppmi");
    get_probs(tp);
    report_pmi(topfile, pmifile, toppmifile, ddN.T, ddN.W, 1, 
               pmicount, tp, tpmi);
    free(toppmifile);
    free(pmifile);
    free(tp);
  }

  /*
   *   now report words and diagnostics
   */
  //ttop_open(topfile);
  if ( fullreport ) {
    rp = fopen(repfile,"w");
    if ( !rp ) 
      yap_sysquit("Cannot open file '%s' for write\n", repfile);
    fprintf(rp, "#topic index rank prop word-sparse doc-sparse eff-words eff-docs docs-bound top-one "
	    "dist-unif dist-unigrm");
    if ( PCTL_BURSTY() ) 
      fprintf(rp, " burst-concent");
    if ( ddN.tokens )  
      fprintf(rp, " ave-length");
    fprintf(rp, " coher");
    if ( pmicount ) 
      fprintf(rp, " pmi");
    fprintf(rp, "\n#word topic index rank");
    if ( ddS.Nwt )
      fprintf(rp, " count");
    fprintf(rp, " prop cumm df coher\n");
    
  }
  for (k=0; k<ddN.T; k++) {
    int cnt, termcnt = 0;
    int kk = psort[k];
    uint32_t **dfmtx;

    if ( ddP.phi==NULL && ddS.NWt[kk]==0 )
      continue;
    /*
     *   grab word prob vec for later use
     */
    if ( ddS.Nwt ) {
      int w;
      for (w=0; w<ddN.W; w++)
	pvec[w] = wordprob(w,kk);
    } else if ( ddP.phi ) 
      fv_copy(pvec, ddP.phi[kk], ddN.W);
    else if ( ddS.phi ) 
      fv_copy(pvec, ddS.phi[kk], ddN.W);

    /*
     *  rebuild word list
     */
    tscorek = kk;
    cnt = buildindk(kk, indk);
    topk(topword, cnt, indk, tscore);
    if ( topword<cnt )
      cnt = topword;
    assert(cnt>0);
    if ( termstats ) {
      termcnt = buildtermindk(kk, termindk, termstats);
      topk(topword, termcnt, termindk, termtscore);
      if ( topword<termcnt )
	termcnt = topword;
    }
    /*
     *     df stats for topic returned as matrix
     */
    dfmtx = hca_dfmtx(indk, cnt, kk);

    if ( ddS.Nwt && (ddS.NWt[kk]*ddN.T*100<Nk_tot || ddS.NWt[kk]<5 )) 
      underused++;
    /*
     *  print stats for topic
     *    Mallet:  tokens, doc_ent, ave-word-len, coher., 
     *             uni-dist, corp-dist, eff-no-words
     */
    yap_message("Topic %d/%d", kk, k);
    {
      /*
       *   compute diagnostics
       */
      double prop = gtvec[kk];
      float *dprop = docprop(kk);
      double spw = 0;
      double spd = ((double)nonzero_Ndt(kk))/((double)ddN.DT);
#ifdef KL
      double ew = fv_kl(p0vec,pvec,ddN.W);
#else
      double ew = exp(fv_entropy(pvec,ddN.W));
#endif
      double ud = fv_helldistunif(pvec,ddN.W);
      double pd = fv_helldist(pvec,gpvec,ddN.W);
      double sl = fv_avestrlen(pvec,ddN.tokens,ddN.W);
      double co = coherence(dfmtx, cnt);
      double ed = dprop?exp(fv_entropy(dprop,ddN.DT)):ddN.DT;
#define MALLET_EW
#ifdef MALLET_EW
      double ewp = dprop?(1.0/fv_expprob(pvec,ddN.W)):ddN.W;
#endif
      double da = dprop?fv_bound(dprop,ddN.DT,1.0/sqrt((double)ddN.T)):0;
      sparsitydoc += spd;
      yap_message((ddN.T>200)?" p=%.3lf%%":" p=%.2lf%%",100*prop);   
      if ( ddS.Nwt ) {
	spw = ((double)nonzero_Nwt(kk))/((double)ddN.W);
	sparsityword += spw;
	yap_message(" ws=%.1lf%%", 100*(1-spw));
      } 
      yap_message(" ds=%.1lf%%", 100*(1-spd) );
#ifdef KL
      yap_message(" ew=%lf", ew);
#else
      yap_message(" ew=%.0lf", ew);
#endif
#ifdef MALLET_EW
      yap_message(" ewp=%.1lf", ewp); 
#endif
      yap_message(" ed=%.1lf", ed); 
      yap_message(" da=%.0lf", da+0.1); 
      yap_message(" t1=%u", top1cnt[kk]); 
      yap_message(" ud=%.3lf", ud); 
      yap_message(" pd=%.3lf", pd); 
      if ( PCTL_BURSTY() ) 
	yap_message(" bd=%.3lf", ddP.bdk[kk]); 
      if ( ddP.NGbeta ) {
	/*
	 *   approx. as sqrt(var(lambda_k)/lambda-normaliser
	 */
	double ngvar = sqrt(ddP.NGalpha[kk])
	  * (ngalpha[kk]/ddP.NGalpha[kk]);
	yap_message(" ng=%.4lf,%.4lf", 
		    ngalpha[kk], ngvar/ngalpha[kk]);
	if ( verbose>2 )
	    yap_message(" ngl=%.4lf,%.4lf, nga=%.4lf,%.4lf", 
		    ddP.NGalpha[kk]/ddP.NGbeta[kk], 
		    sqrt(ddP.NGalpha[kk]/ddP.NGbeta[kk]/ddP.NGbeta[kk]),
		    ddP.NGalpha[kk], ddP.NGbeta[kk]); 
      }
      if ( ddN.tokens )  
	yap_message(" sl=%.2lf", sl); 
      yap_message(" co=%.3lf%%", co);
      if ( pmicount ) 
	yap_message(" pmi=%.3f", tpmi[kk]);
      if ( fullreport ) {
	fprintf(rp,"topic %d %d", kk, k);
	fprintf(rp," %.6lf", prop);   
	if ( ddS.Nwt ) {
	  fprintf(rp," %.6lf", (1-spw));
	} else {
	  fprintf(rp," 0");
	}
	fprintf(rp," %.6lf", (1-spd) );
#ifdef KL
	yap_message(" %lf", ew);
#else
	fprintf(rp," %.2lf", ew);
#endif
#ifdef MALLET_EW
	fprintf(rp," %.2lf", ewp); 
#endif
	fprintf(rp," %.2lf", ed); 
	fprintf(rp," %.0lf", da+0.1); 
	fprintf(rp," %u", top1cnt[kk]); 
	fprintf(rp," %.6lf", ud); 
	fprintf(rp," %.6lf", pd); 
	if ( PCTL_BURSTY() ) 
	  fprintf(rp," %.3lf", ddP.bdk[kk]); 
	fprintf(rp," %.4lf", (ddN.tokens)?sl:0); 
	fprintf(rp," %.6lf", co);
	if ( pmicount ) 
	  fprintf(rp," %.4f", tpmi[kk]);
	fprintf(rp,"\n");
      }
      if ( dprop) free(dprop);
    }
    if ( verbose>1 ) {
      double pcumm = 0;
      /*
       *   print top words:
       *     Mallet:   rank, count, prob, cumm, docs, coh
       */
      yap_message("\ntopic %d/%d", kk, k);
      yap_message(" words=");
      for (w=0; w<cnt; w++) {
	if ( w>0 ) yap_message(",");
	if ( ddN.tokens ) 
	  yap_message("%s", ddN.tokens[indk[w]]);
	else
	  yap_message("%d", indk[w]);
	if ( verbose>2 ) {
	  if ( scoretype == ST_count )
	    yap_message("(%d)", (int)(tscore(indk[w])+0.2));
	  else
	    yap_message("(%6lf)", tscore(indk[w]));
	}
	if ( fullreport ) {
	  fprintf(rp, "word %d %d %d", kk, indk[w], w);
	  if ( ddS.Nwt )
	    fprintf(rp, " %d", ddS.Nwt[indk[w]][kk]);
	  pcumm += pvec[indk[w]];
	  fprintf(rp, " %.6f %.6f", pvec[indk[w]], pcumm);
	  fprintf(rp, " %d", dfmtx[w][w]); 
	  fprintf(rp, " %.6f", coherence_word(dfmtx, cnt, w));
	  if ( ddN.tokens ) 
	    fprintf(rp, " %s", ddN.tokens[indk[w]]);
	  fprintf(rp, "\n");
	}
      }
      if ( termstats ) {
	yap_message(" terms=");
	for (w=0; w<termcnt; w++) {
	  if ( w>0 ) yap_message(",");
	  if ( ddN.tokens ) 
	    yap_message("%s", termstats->tokens[termindk[w]]);
	  else
	    yap_message("%d", termstats->Kmin+termindk[w]);
	  if ( verbose>2 ) {
	    if ( scoretype == ST_count )
	      yap_message("(%d)", (int)(termtscore(termindk[w])+0.2));
	    else
	      yap_message("(%6lf)", termtscore(termindk[w]));
	  }
	  if ( fullreport ) {
	    fprintf(rp, "term %d %d %d", kk, termindk[w], w);
	    fprintf(rp, " %d", termstats->Nkt[termindk[w]][kk]);
	    fprintf(rp, " %s", termstats->tokens[termindk[w]]);
	    fprintf(rp, "\n");
	  }
	}
      }
    }
    yap_message("\n");
    free(dfmtx[0]); free(dfmtx); 
  }
  if ( verbose>1 && ddP.PYbeta && (ddP.phi==NULL || ddP.betapr) ) {
    int cnt;
    double pcumm = 0;
     /*
     *    print root words
     */
    tscorek = -1;
    cnt = buildindk(-1,indk);
    topk(topword, cnt, indk, (ddP.phi==NULL)?countscore:phiscore);
    /*
     *     cannot build df mtx for root because
     *     it is latent w.r.t. topics
     */
    yap_message("Topic root words=");
    if ( fullreport ) {
      int w;
      for (w=0; w<ddN.W; w++)
	pvec[w] = betabasewordprob(w);
#ifdef KL
      double ew = fv_kl(p0vec,pvec,ddN.W);
#else
      double ew = exp(fv_entropy(pvec,ddN.W));
#endif

      double ud = fv_helldistunif(pvec,ddN.W);
      double pd = fv_helldist(pvec,gpvec,ddN.W);
      fprintf(rp,"topic -1 -1 0 0");
      fprintf(rp," %.4lf", ew); 
      fprintf(rp," %.6lf", ud); 
      fprintf(rp," %.6lf", pd); 
      fprintf(rp,"\n");
    }
    for (w=0; w<topword && w<cnt; w++) {
      if ( w>0 ) yap_message(",");
      if ( ddN.tokens )
	yap_message("%s", ddN.tokens[indk[w]]);
      else
	yap_message("%d", indk[w]);
      if ( verbose>2 )
	yap_message("(%6lf)", countscore(indk[w]));
      if ( fullreport ) {
	fprintf(rp, "word %d %d %d", -1, indk[w], w);
	if ( ddS.TwT )
	  fprintf(rp, " %d", ddS.TwT[w]);
	pcumm += pvec[indk[w]];
	fprintf(rp, " %.6f %.6f", pvec[indk[w]], pcumm);
	fprintf(rp, " 0 0"); 
	if ( ddN.tokens ) 
	  fprintf(rp, " %s", ddN.tokens[indk[w]]);
	fprintf(rp, "\n");
      }   
    }
    yap_message("\n");
  }
  yap_message("\n");
  if ( rp )
    fclose(rp);
	     
  if ( ddS.Nwt )
    yap_message("Average topicXword sparsity = %.2lf%%\n",
                100*(1-sparsityword/ddN.T) );
  yap_message("Average docXtopic sparsity = %.2lf%%\n"
	      "Underused topics = %.1lf%%\n",
	      100*(1-sparsitydoc/ddN.T), 
	      100.0*underused/(double)ddN.T);
  if ( pmicount ) 
    yap_message("Average PMI = %.3f\n", tpmi[ddN.T]);

  /*
   *   print 
   */
  if ( 1 ) {
    float **cmtx = hca_topmtx();
    int t1, t2;
    int m1, m2;
    float mval;
    char *corfile = yap_makename(resstem,".topcor");
    fp = fopen(corfile,"w");
    if ( !fp ) 
      yap_sysquit("Cannot open file '%s' for write\n", corfile);
   /*
    *   print file
     */
    for (t1=0; t1<ddN.T; t1++) {
      for (t2=0; t2<t1; t2++) 
	 if ( cmtx[t1][t2]>1.0e-7 ) 
	  fprintf(fp, "%d %d %0.6f\n", t1, t2, cmtx[t1][t2]);
    }
    fclose(fp);
    free(corfile);
    /*
     *   display maximum
     */
    m1 = 1; m2 = 0;
    mval = cmtx[1][0];
    for (t1=0; t1<ddN.T; t1++) {
      for (t2=0; t2<t1; t2++) {
	if ( mval<cmtx[t1][t2] ) {
	  mval = cmtx[t1][t2];
	  m1 = t1;
	  m2 = t2;
	}
      }
    }
    yap_message("Maximum correlated topics (%d,%d) = %f\n", m1, m2, mval);
    free(cmtx[0]); free(cmtx);
  }

  /*
   *  print burstiness report
   */
  if ( PCTL_BURSTY() ) {
    int tottbl = 0;
    int totmlttbl = 0;
    int totmlt = 0;
    int i;
    for (i=0; i<ddN.NT; i++) {
      if ( Z_issetr(ddS.z[i]) ) {
	if ( M_multi(i) )
	  totmlttbl++;
	tottbl++;
      }
      if ( M_multi(i) )
	totmlt++;
    }
    yap_message("Burst report: multis=%.2lf%%, tables=%.2lf%%, tbls-in-multis=%.2lf%%\n",
		100.0*((double)ddM.dim_multiind)/ddN.N,
		100.0*((double)tottbl)/ddN.NT,
		100.0*((double)totmlttbl)/totmlt);
  }
  yap_message("\n");

  free(topfile);
  if ( repfile ) free(repfile);
  if ( top1cnt ) free(top1cnt);
  free(indk);
  free(psort);
  if ( ngalpha )
    free(ngalpha);
  if ( pmicount )
    free(tpmi);
  if ( NwK ) {
    free(NwK);
    NwK = NULL;
  }
#ifdef KL
  free(ddD.df);
#endif
  free(pvec); 
  free(gtvec);
  free(gpvec);
  tstats_free(termstats);
}
