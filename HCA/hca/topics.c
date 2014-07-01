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

/*************************************************
 *  return sort order for topics by size
 */
static int pcompar(const void *a, const void *b) {
  uint32_t na = ddS.NWt[*(uint32_t*)a];
  uint32_t nb = ddS.NWt[*(uint32_t*)b];
  if ( na<nb ) return 1;
  if ( na>nb ) return -1;
  return 0;
}
static uint32_t *sorttops(int K) {
  uint32_t *psort = u32vec(K);
  int k;
  for (k=0; k<K; k++) psort[k] = k;
  qsort(psort, K, sizeof(*psort), pcompar);
  return psort;
}

/*************************************************
 *  find top-K by bubble sort
 */
static void topk(int K, int N, int *ind, double (*foo)(int)) {
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


/***************************************************
 *
 *   scoring utilities
 */
/*
 *  for IDF score
 */
static uint32_t NWK = 0;
static uint32_t *NwK = NULL;

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

static double idfscore(int w) {
  return (getn(w)+0.2)/(NwK[w]+0.2*ddN.T);
}

static double phiscore(int w) {
  assert(ddP.phi || ddS.phi);
  if ( tscorek<0 ) 
    return ddS.TwT[w];
  if ( ddP.phi )
    return ddP.phi[tscorek][w];
  return ddS.phi[tscorek][w];
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

/********************************************************
 *  build non-zero indices for topic
 */
static int buildindk(int k, int *indk) {
  int w;
  int cnt=0;
  if ( k<0 ) {
    /*
     *  for root topic
     */
    if ( ddP.phi!=NULL ) 
      return 0;
    for (w=0; w<ddN.W; w++) {
      if ( ddS.TwT[w]>0 ) indk[cnt++] = w;
    }
    return cnt;
  }
  if ( ddP.phi==NULL ) {
    for (w=0; w<ddN.W; w++) {
      if ( ddS.Nwt[w][k]>0 ) indk[cnt++] = w;
    }
  } else {
    float **phi;
    if ( ddP.phi )
      phi = ddP.phi;
    else
      phi = ddS.phi;
    for (w=0; w<ddN.W; w++) {
      if ( phi[k][w]>0.5/ddN.W ) indk[cnt++] = w;
    }
  }
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
                       enum ScoreType scoretype, int pmicount) {
  int w,k;
  int *indk = NULL;
  int Nk_tot = 0;
  double (*tscore)(int) = NULL;
  double sparsityword = 0;
  double sparsitydoc = 0;
  double underused = 0;
  int nophi = (ddP.phi==NULL) && (ddS.phi==NULL);
  FILE *fp;
  float *tpmi;
  char *topfile;
  uint32_t *psort;

  /*
   *  topic stats
   */
  double *spw;
  double *spd; 
  float *effcnt;

  if ( pmicount>topword )
    pmicount = topword;
  if ( scoretype == ST_idf ) {
    tscore = idfscore;
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

  /*
   *  first collect counts of each word/term
   */
  if ( scoretype != ST_count && scoretype != ST_phi ) {
    NwK = u32vec(ddN.W);
    if ( !NwK )
      yap_quit("Out of memory in hca_displaytopics()\n");
    for (w=0; w<ddN.W; w++) {
      NwK[w] = 0;
    }
    NWK = 0;
    for (w=0; w<ddN.W; w++) {
      for (k=0; k<ddN.T; k++) {
	NwK[w] += ddS.Nwt[w][k];    //  should use CCT_ReadN()
      }
      NWK += NwK[w];
    }
  }

  for (k=0; k<ddN.T; k++) {
    Nk_tot += ddS.NWt[k];
  }

  psort = sorttops(ddN.T);

  if ( pmicount ) {
    tpmi = malloc(sizeof(*tpmi)*(ddN.T+1));
    if ( !tpmi )
      yap_quit("Cannot allocate tpmi in hca_displaytopics()\n");
  }
  spd = malloc(sizeof(*spd)*ddN.T);
  spw = malloc(sizeof(*spw)*ddN.T);
  if ( !spw || !spd )
    yap_quit("Cannot allocate spw,spd in hca_displaytopics()\n");
  if ( !nophi ) {
    effcnt = malloc(sizeof(*effcnt)*ddN.T);
    if ( !effcnt )
      yap_quit("Cannot allocate effcnt in hca_displaytopics()\n");
  }
  indk = malloc(sizeof(*indk)*ddN.W);
  if ( !indk )
    yap_quit("Cannot allocate indk in hca_displaytopics()\n");

  /*
   *   two passes through, 
   *           first to build the top words and dump to file
   */
  topfile = yap_makename(resstem,".top");
  fp = fopen(topfile,"w");
  if ( !fp ) 
    yap_sysquit("Cannot open file '%s' for write\n", topfile);
  yap_message("\n");
  for (k=0; k<ddN.T; k++) {
    int cnt;
    tscorek = k;
    /*
     *    build sorted word list
     */
    cnt = buildindk(k, indk);
    topk(topword, cnt, indk, tscore);
    /*
     *   compute diagnostics
     */
    spd[k] = ((double)nonzero_Ndt(k))/((double)ddN.DT);
    sparsitydoc += spd[k];
    if ( nophi ) {
      spw[k] = ((double)nonzero_Nwt(k))/((double)ddN.W);
      sparsityword += spw[k];
    }
    if ( ddS.NWt[k]*ddN.T*100<Nk_tot || ddS.NWt[k]<5 ) 
      underused++;
    if ( !nophi )
      effcnt[k] = exp(phi_entropy(k));
    if ( cnt==0 )
      continue;
    /*
     *   dump words to file
     */
    fprintf(fp,"%d: ", k);
    for (w=0; w<topword && w<cnt; w++) {
      fprintf(fp," %d", (int)indk[w]);
    }
    fprintf(fp, "\n");
  }
  if ( ddP.PYbeta && nophi ) {
    int cnt;
     /*
     *    dump root words
     */
    tscorek = -1;
    cnt = buildindk(-1, indk);
    topk(topword, cnt, indk, countscore);
    fprintf(fp,"-1:");
    for (w=0; w<topword && w<cnt; w++) {
      fprintf(fp," %d", (int)indk[w]);
    }
    fprintf(fp, "\n");
  }
  fclose(fp);
  if ( verbose>1 ) yap_message("\n");

  if ( pmicount ) {
    char *toppmifile;
    char *pmifile;
    double *tp;
    tp = dvec(ddN.T);
    pmifile=yap_makename(stem,".pmi");
    toppmifile=yap_makename(topfile,"pmi");
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
  ttop_open(topfile);
  for (k=0; k<ddN.T; k++) {
    int cnt;
    int kk = psort[k];
    if ( ddS.NWt[kk]==0 )
      continue;
    /*
     *  rebuild word list
     */
    tscorek = kk;
    cnt = buildindk(kk, indk);
    topk(topword, cnt, indk, tscore);
    /*
     *  print stats
     */
    yap_message("Topic %d/%d", kk, k);
    if ( ddP.phi==NULL ) 
      yap_message((ddN.T>200)?" p=%.3lf%%,":" p=%.2lf%%,", 
		  100*((double)ddS.NWt[kk])/(double)Nk_tot);   
    if ( nophi ) 
      yap_message(" ws=%.1lf%%,", 100*(1-spw[kk]));
    else
      yap_message(" #=%.0lf,", effcnt[kk]); 
    yap_message(" ds=%.1lf%%", 100*(1-spd[kk]) );
    if ( pmicount ) 
      yap_message(" pmi=%.3f,", tpmi[kk]);
    if ( verbose>1 ) {
      /*
       *   print top words
       */
      yap_message(" words=");
      for (w=0; w<topword && w<cnt; w++) {
	if ( w>0 ) yap_message(",");
	if ( ddN.tokens ) 
	  yap_message("%s", ddN.tokens[indk[w]]);
	else
	  yap_message("%d", indk[w]);
	if ( verbose>2 )
	  yap_message("(%6lf)", tscore(indk[w]));
      }
    }
    yap_message("\n");
  }
  if ( verbose>1 && ddP.PYbeta && nophi ) {
    int cnt;
     /*
     *    print root words
     */
    tscorek = -1;
    cnt = buildindk(-1,indk);
    topk(topword, cnt, indk, countscore);
    yap_message("Topic root words =");
    for (w=0; w<topword && w<cnt; w++) {
      if ( w>0 ) yap_message(",");
      if ( ddN.tokens )
	yap_message("%s", ddN.tokens[indk[w]]);
      else
	yap_message("%d", indk[w]);
      if ( verbose>2 )
	yap_message("(%6lf)", countscore(indk[w]));
    }
    yap_message("\n");
  }
  yap_message("\n");
	     
  if ( nophi )
    yap_message("Average topicXword sparsity = %.2lf%%\n",
                100*(1-sparsityword/ddN.T) );
  yap_message("Average docXtopic sparsity = %.2lf%%\n"
	      "Underused topics = %.1lf%%\n",
	      100*(1-sparsitydoc/ddN.T), 
	      100.0*underused/(double)ddN.T);
  if ( pmicount ) 
    yap_message("Average PMI = %.3f\n", tpmi[ddN.T]);

  /*
   *  print burstiness report
   */
  if ( ddP.bdk!=NULL) {
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
  free(indk);
  free(spd);
  free(spw);
  free(psort);
  if ( !nophi )
    free(effcnt);
  if ( pmicount )
    free(tpmi);
  if ( scoretype != ST_count ) {
    free(NwK);
    NwK = NULL;
  }
}
