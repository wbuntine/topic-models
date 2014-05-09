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

/*
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


/*
 *  for IDF score
 */
static uint32_t NWK = 0;
static uint32_t *NwK = NULL;

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
  double N = NwK[w];
  return getn(w)/((double)ddS.NWt[tscorek]) 
    - (N*ddS.NWt[tscorek]/((double)ddN.NT))/((double)ddN.NT);
}

static double countscore(int w) {
  return getn(w);
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

void hca_displaytopics(char *resstem, int topword, enum ScoreType scoretype) {
  int w,k;
  int *indk = NULL;
  int Nk_tot = 0;
  double (*tscore)(int) = NULL;
  double sparsityword = 0;
  double sparsitydoc = 0;
  double underused = 0;
  char *fname = yap_makename(resstem,".top");
  int nophi = (ddP.phi==NULL) && (ddS.phi==NULL);
  FILE *fp;

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

  fp = fopen(fname,"w");
  if ( !fp ) 
    yap_sysquit("Cannot open file '%s' for write\n", fname);

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

  assert(ddN.tokens);

  for (k=0; k<ddN.T; k++) {
    Nk_tot += ddS.NWt[k];
  }

  indk = malloc(sizeof(*indk)*ddN.W);
  if ( !indk )
    yap_quit("Cannot allocate indk\n");
  
  for (k=0; k<ddN.T; k++) {
    int cnt;
    double spw;
    double spd; 
    tscorek = k;
    /*
     *    print top words
     */
    cnt=0;
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
    topk(topword, cnt, indk, tscore);
    spd = ((double)nonzero_Ndt(k))/((double)ddN.DT);
    sparsitydoc += spd;
    if ( nophi ) {
      spw = ((double)nonzero_Nwt(k))/((double)ddN.W);
      sparsityword += spw;
    }
    if ( ddS.NWt[k]*ddN.T*100<Nk_tot ) 
      underused++;
    yap_message("\nTopic %d (", k);
    if ( ddP.phi==NULL ) 
      yap_message("p=%.2lf%%,", 
		  100*((double)ddS.NWt[k])/(double)Nk_tot);   
    if ( nophi ) 
      yap_message("ws=%.1lf%%,", 100*(1-spw));
    else
      yap_message("#=%.0lf,", exp(phi_entropy(k)));
    yap_message("ds=%.1lf%%", 100*(1-spd) );
    fprintf(fp,"%d: ", k);
    yap_message(") words =");
    for (w=0; w<topword && w<cnt; w++) {
      fprintf(fp," %d", (int)indk[w]);
      if ( verbose>2 ) {
	double score = tscore(indk[w]);
	yap_message(",%s(%6lf)", ddN.tokens[indk[w]], score);
      } else
	yap_message(",%s", ddN.tokens[indk[w]]);
    }
    yap_message("\n");
    fprintf(fp, "\n");
  }
	     
  if ( ddP.PYbeta && nophi ) {
    int cnt;
     /*
     *    print root words
     */
    tscorek = -1;
    cnt=0;
    for (w=0; w<ddN.W; w++) {
      if ( ddS.TwT[w]>0 ) indk[cnt++] = w;
    }
    topk(topword, cnt, indk, tscore);
    yap_message("\nTopic root words =");
    fprintf(fp,"-1:");
    for (w=0; w<topword && w<cnt; w++) {
      fprintf(fp," %d", (int)indk[w]);
      if ( verbose>2 ) {
	double score = tscore(indk[w]);
	yap_message(",%s(%6lf)", ddN.tokens[indk[w]], score);
      } else
	yap_message(",%s", ddN.tokens[indk[w]]);
    }
    yap_message("\n");
    fprintf(fp, "\n");
  }
  if ( nophi )
    yap_message("Average topicXword sparsity = %.2lf%%, ",
		100*(1-sparsityword/ddN.T) );
  yap_message("Average docXtopic sparsity = %.2lf%%, "
	      "underused topics = %.1lf%%\n",
	      100*(1-sparsitydoc/ddN.T), 
	      100.0*underused/(double)ddN.T);
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
    yap_message("doc PYP report:   multis=%.2lf%%,  tables=%.2lf%%, tbls-in-multis=%.2lf%%\n",
		100.0*((double)ddM.dim_multiind)/ddN.N,
		100.0*((double)tottbl)/ddN.NT,
		100.0*((double)totmlttbl)/totmlt);
  }
  fclose(fp);
  free(fname);
  free(indk);
  if ( scoretype != ST_count ) {
    free(NwK);
    NwK = NULL;
  }
}
