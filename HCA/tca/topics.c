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
 *  for IDF score and elsewhere, these are totals built using
 *   getnk() for each epoch
 */
static uint32_t NWK = 0;
static uint32_t *NwK = NULL;

/*
 *    when sorting, need to know these indices,
 *    so we make them statics, ...
 */
static int tscorek;
static int tscoree;

static unsigned getnk(int w, int k) {
  if ( ddP.phi && ddP.mu )
    return  round(ddP.mu[tscoree][k] * ddP.phi[tscoree][w][k] * NWK);
  return ddS.m_vte[w][k][tscoree]
    + ((ddP.phi==NULL&&tscoree+1<ddN.E)?ddS.s_vte[w][k][tscoree+1]:0);
}
static unsigned getNk(int k) {
  if ( ddP.phi && ddP.mu )
    return round(ddP.mu[tscoree][k] * NWK);
  return ddS.M_Vte[k][tscoree]
    + ((ddP.phi==NULL&&tscoree+1<ddN.E)?ddS.S_Vte[k][tscoree+1]:0);
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
  int M_tot = 0;
  double (*tscore)(int) = NULL;
  double sparsityword = 0;
  double sparsitydoc = 0;
  double underused = 0;
  char *fname = yap_makename(resstem,".top");
  FILE *fp;
  /*
   *  store mu and phi estimates computed iteratively
   */
  double *pvec;
  double **wmtx = NULL;

  if ( ! ddP.mu )
    pvec = dvec(ddN.T);
  if ( !ddP.phi ) 
    wmtx = dmat(ddN.W, ddN.T);

  assert(ddN.tokens);
  if ( !ddP.phi && ! wmtx )
    yap_quit("Cannot allocate memory tca_displaytopics()\n"); 

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

  indk = malloc(sizeof(*indk)*ddN.W);
  if ( !indk )
    yap_quit("Cannot allocate memory tca_displaytopics()\n"); 
  NwK = u32vec(ddN.W);
  if ( !NwK )
    yap_quit("Out of memory in cca_displaytopics()\n");
  fp = fopen(fname,"w");
  if ( !fp ) 
    yap_sysquit("Cannot open file '%s' for write\n", fname);

  /*
   *   initialise pvec and wmtx since computed iteratively
   */
  if ( ! ddP.mu )
    mu_prob_iter(-1, pvec);
  if ( ! ddP.phi )
    phi_prob_iter(-1, wmtx);

  for (tscoree=0; tscoree<ddN.E; tscoree++) {
    /*
     *   update estimates for mu and phi
     */
    if ( ! ddP.mu )
      mu_prob_iter(tscoree, pvec);
    if ( ! ddP.phi )
      phi_prob_iter(tscoree, wmtx);

    /*
     *  first collect counts of each word/term in epoch
     */
    if ( !ddP.phi ) {
      /*
       *   get from stats
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
	M_tot += ddS.M_Vte[k][tscoree];
      }
    } else {
      /*
       *  get from raw data
       */
      int e, l;
      int endd, startd = 0;
      for (w=0; w<ddN.W; w++)
	NwK[w] = 0;
      for (e=0; e<tscoree; e++)
	startd += ddD.esize[e];
      endd = startd + ddD.esize[tscoree];
      for (l=ddD.N_dTcum[startd]; l<ddD.N_dTcum[endd]; l++)
	NwK[ddD.w[l]]++;
      NWK = 0;
      for (w=0; w<ddN.W; w++)
	NWK += NwK[w];
    }
    
    yap_message("\nEpoch %d\n", (int)tscoree);
    for (k=0; k<ddN.T; k++) {
      int cnt;
      double spw=0;
      double spd=0; 
      tscorek = k;
      /*
       *    print top words
       */
      cnt=0;
      for (w=0; w<ddN.W; w++) {
	if ( ddS.m_vte[w][k][tscoree]>0 ) indk[cnt++] = w;
      }
      topk(topword, cnt, indk, tscore);
      if ( !ddP.phi ) {
	spd = 1-((double)nonzero_n_dt(k))/((double)ddN.DT);
	sparsitydoc += spd;
	spw = 1-((double)nonzero_m_vte(tscoree,k))/((double)ddN.W);
	sparsityword += spw;
	if ( !ddP.phi && ddS.M_Vte[k][tscoree]*ddN.T*100<M_tot ) 
	  underused++;
      }

      yap_message("Topic %d/%d (", k, (int)tscoree);
      if ( !ddP.mu ) 
	yap_message("p=%.2lf%%/%.2lf%%", 
		    100.0*((double)ddS.M_Vte[k][tscoree])/(double)M_tot,
		    100.0*pvec[k]);
      else
	yap_message("p=%.2lf%%", 100.0*ddP.mu[tscoree][k]);
      if ( !ddP.phi ) {      
	yap_message(",ws=%.1lf%%,", 100*spw);
	yap_message("ds=%.1lf%%", 100*spd );
      }
      fprintf(fp,"%d,%d: ", (int)tscoree, (int)k);
      yap_message(") words =");
      for (w=0; w<topword && w<cnt; w++) {
        double prob;
        if ( verbose>3 )
          prob = ddS.m_vte[indk[w]][k][tscoree]/(double)ddS.M_Vte[k][tscoree];
	/*  print to file */
	fprintf(fp," %d", (int)indk[w]);
	if ( verbose>2 ) {
	  fprintf(fp,"(");
	  fprintf(fp,(scoretype == ST_count)?"%.0lf":"%6lf",tscore(indk[w]));
	  if ( verbose>3 ) 
	    fprintf(fp,",%6lf/%6lf", prob, 
                    ddP.phi?ddP.phi[tscoree][indk[w]][k]:wmtx[indk[w]][k]);
	  fprintf(fp,")");
	} 	
	/*  yap to report */
	yap_message(",%s", ddN.tokens[indk[w]]);
	if ( verbose>2 ) {
	  yap_message("(");
	  if ( verbose>3 ) yap_message("s=");
	  yap_message((scoretype == ST_count)?"%.0lf":"%6lf",tscore(indk[w]));
	  if ( verbose>3 )
	    yap_message(",p=%6lf/%6lf", prob, 
                        ddP.phi?ddP.phi[tscoree][indk[w]][k]:wmtx[indk[w]][k]);
	  yap_message(")");
	} 
      }
      yap_message("\n");
      fprintf(fp, "\n");
    }
  }

  if ( !ddP.phi ) {
    yap_message("Average topicXword sparsity = %.2lf%%, ",
		100*(sparsityword/ddN.T/ddN.E) );
    yap_message("Average docXtopic sparsity = %.2lf%%, "
		"underused topics = %.1lf%%\n",
		100*(sparsitydoc/ddN.T/ddN.E), 
		100.0*underused/ddN.T/ddN.E);
  }
  
  fclose(fp);
  if ( !ddP.mu ) 
    free(pvec);
  if ( !ddP.phi ) {
    free(wmtx[0]);  free(wmtx);
  }
  free(fname);
  free(indk);
  free(NwK);
  NwK = NULL;
}
