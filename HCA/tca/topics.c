/*
 * Topic reporting routines
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
#include "tca.h"
#include "stats.h"
#include "lgamma.h" 
#include "util.h" 
#include "yap.h" 

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
static int tscoree;

static unsigned getnk(int w, int k) {
  return ddS.m_evt[tscoree][w][k]
    + ((ddP.phi==NULL&&tscoree+1<ddN.E)?ddS.s_evt[tscoree+1][w][k]:0);
}
static unsigned getNk(int k) {
  return ddS.M_eVt[tscoree][k]
    + ((ddP.phi==NULL&&tscoree+1<ddN.E)?ddS.S_eVt[tscoree+1][k]:0);
}
static unsigned getn(int w) {
  return getnk(w,tscorek);
}
static unsigned getN() {
  return getNk(tscorek);
}

static double idfscore(int w) {
  return (getn(w)+0.2)/(NwK[w]+0.2*ddN.T);
}

static double lowerQ;
static double Qscore(int w) {
  int n = getn(w);
  double N;
  N = NwK[w];
  if ( n/((double)N) <= lowerQ )
    return 0;
  return N/NWK * (n/((double)N)-lowerQ)/(1-lowerQ);
}

static double costscore(int w) {
  double Nk = getN();
  return getn(w)/Nk - (NwK[w]*Nk/(double)NWK/(double)NWK);
}

static double countscore(int w) {
  return getn(w);
}

void tca_displaytopics(char *resstem, int topword, enum ScoreType scoretype) {
  int w,k;
  int *indk = NULL;
  int M_tot;
  double (*tscore)(int) = NULL;
  double sparsityword = 0;
  double sparsitydoc = 0;
  double underused = 0;
  char *fname = yap_makename(resstem,".top");
  FILE *fp;

  assert(ddN.tokens);
  
  if ( scoretype == ST_idf ) {
    tscore = idfscore;
  } else if ( scoretype == ST_count ) {
    tscore = countscore;
  } else if ( scoretype == ST_cost ) {
    tscore = costscore;
  } else if ( scoretype == ST_Q ) {
    tscore = Qscore;
    lowerQ = 1.0/ddN.T;
  }    

  NwK = u32vec(ddN.W);
  if ( !NwK )
    yap_quit("Out of memory in cca_displaytopics()\n");
  fp = fopen(fname,"w");
  if ( !fp ) 
    yap_sysquit("Cannot open file '%s' for write\n", fname);

  for (tscoree=0; tscoree<ddN.E; tscoree++) {
    /*
     *  first collect counts of each word/term
     */
    NWK = 0;
    for (w=0; w<ddN.W; w++) {
      NwK[w] = 0;
      for (k=0; k<ddN.T; k++) {
	NwK[w] += getnk(w,k); 
      }
      NWK += NwK[w];
    }    
    M_tot = 0;
    for (k=0; k<ddN.T; k++) {
      M_tot += ddS.M_eVt[tscoree][k];
    }
    
    indk = malloc(sizeof(*indk)*ddN.W);
    if ( !indk )
      yap_quit("Cannot allocate indk\n");
    
    yap_message("\nEpoch %d\n", (int)tscoree);
    for (k=0; k<ddN.T; k++) {
      int cnt;
      double spw;
      double spd; 
      tscorek = k;
      /*
       *    print top words
       */
      cnt=0;
      for (w=0; w<ddN.W; w++) {
	if ( ddS.m_evt[tscoree][w][k]>0 ) indk[cnt++] = w;
      }
      topk(topword, cnt, indk, tscore);
      spd = ((double)nonzero_n_dt(k))/((double)ddN.DT);
      if ( tscoree==0 )
	sparsitydoc += spd;
      spw = ((double)nonzero_m_evt(tscoree,k))/((double)ddN.W);
      sparsityword += spw;
      
      if ( ddS.M_eVt[tscoree][k]*ddN.T*100<M_tot ) 
	underused++;
      yap_message("Topic %d/%d (", k, (int)tscoree);
      yap_message("p=%.2lf%%/%.2lf%%,", 
		  100.0*((double)ddS.M_eVt[tscoree][k])/(double)M_tot,
		  100.0*getNk(k)/(double)NWK);
      yap_message("ws=%.1lf%%,", 100*(1-spw));
      yap_message("ds=%.1lf%%", 100*(1-spd) );
      fprintf(fp,"%d,%d: ", (int)tscoree, (int)k);
      yap_message(") words =");
      for (w=0; w<topword && w<cnt; w++) {
	/*  print to file */
	fprintf(fp," %d", (int)indk[w]);
	if ( verbose>2 ) {
	  fprintf(fp,"(");
	  fprintf(fp,(scoretype == ST_count)?"%.0lf":"%6lf",tscore(indk[w]));
	  if ( verbose>3 ) {
	    double prob = 
	      ddS.m_evt[tscoree][indk[w]][k]/(double)ddS.M_eVt[tscoree][k];
	    fprintf(fp,",%6lf", prob);
	  } 
	  fprintf(fp,")");
	} 	
	/*  yap to report */
	yap_message(",%s", ddN.tokens[indk[w]]);
	if ( verbose>2 ) {
	  yap_message("(");
	  if ( verbose>3 ) yap_message("s=");
	  yap_message((scoretype == ST_count)?"%.0lf":"%6lf",tscore(indk[w]));
	  if ( verbose>3 ) {
	    double prob = 
	      ddS.m_evt[tscoree][indk[w]][k]/(double)ddS.M_eVt[tscoree][k];
	    yap_message(",p=%6lf", prob);
	  } 
	  yap_message(")");
	} 
      }
      yap_message("\n");
      fprintf(fp, "\n");
    }
  }

  yap_message("Average topicXword sparsity = %.2lf%%, ",
	      100*(1-sparsityword/ddN.T) );
  yap_message("Average docXtopic sparsity = %.2lf%%, "
	      "underused topics = %.1lf%%\n",
	      100*(1-sparsitydoc/ddN.T), 
	      100.0*underused/(double)ddN.T);
  fclose(fp);
  free(fname);
  free(indk);
  free(NwK);
  NwK = NULL;
}
