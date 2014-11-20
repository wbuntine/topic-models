/**
 * Main driver
 * Copyright (C) 2011-2014 Wray Buntine
 *           (C) 2014 Swapnil Mishra
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@monash.edu)
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <limits.h>

#include "yap.h"
#include "util.h"
#include "stable.h"
#include "lgamma.h"
#include "hca.h"
#include "pctl.h"
#include "data.h"
#include "misi.h"
#include "stats.h"
#include "sample.h"
#include "probs.h"
#include "diag.h"
#include "check.h"

/*
 *  not always threading, but code uses this anyway
 */
#include "pargs.h"

#ifdef H_THREADS
#include <pthread.h>
#endif
#include "atomic.h"

void hca_displaytopics(char *stem, char *resstem, int topword, 
		       enum ScoreType score, int pmicount, int fullreport);
void hca_displayclass(char *resstem);

//==================================================
// global variables
//==================================================

rngp_t rngp = NULL;
int verbose = 0;

// #define QUERY

/*
 *    Dimensions
 */
D_dims_t ddN;
/*
 *  hyperparameters
 */
D_pars_t ddP;
D_pctl_t ddT[ParBeta+1];

/*
 *    Basic data about documents:
 *       we assume maximum number of words in corpus
 *       fits in 32 bit
 */
D_data_t ddD;
/*
 *    Statistics
 */
D_stats_t ddS;
/*
 *    cache
 */
D_cache_t ddC;
/*
 *    diagnostics
 */
D_diag_t ddG;
/*
 *   bursty data stuctures
 */
D_DMi_t ddM;

/*
 *    describes the sampling type
 */
static const char *stype() {
  if ( ddP.PYalpha ) {
    if ( ddP.PYbeta )
      return 
	"H.Pitman-Yor sampler for topics"
	", H.Pitman-Yor sampler for words"
	"\n";
      else
	return 
	  "H.Pitman-Yor sampler for topics"
	  ", Dirichlet sampler for words"
	  "\n";
  } else {
    if ( ddP.PYbeta )
      return 
	"Dirichlet sampler for topics"
	", H.Pitman-Yor sampler for words"
	"\n";
    else
      return 
	"Dirichlet sampler for topics"
	", Dirichlet sampler for words"
	"\n";
  }
}

static void usage() {
  fprintf(stderr,"Commandline:  OPTION+ STEM RESSTEM\n"
	  "  (reads STEM.dit and STEM.wit, results to RESSTEM.*)\n  ");
  fprintf(stderr," Version " HCA_VERSION ", "); 
#ifdef H_THREADS
  fprintf(stderr,"threads, ");
#endif
  fprintf(stderr, "%s", stype());
  fprintf(stderr,
          "  OPTION is choice of:\n"
	  "  setting hyperparameters:\n"
          "   -A/B val[,file] #  alpha/beta Dir. prior and mean is val \n"
          "   -A/B dir[,file] #  alpha/beta Dir. prior and mean is default \n"
          "   -A/B hdp[,file] #  for alpha/beta prior use DP with mean 'file'\n"
          "   -A/B pdp[,file] #  for alpha/beta prior use mean 'file'\n"
	  "      file = 'uniform' (default), 'data' or a filename\n"
          "   -A/B hpdd      #  for alpha/beta prior use truncated GEM\n"
          "   -S var=value   #  initialise var=a,b,a0,b0,aw,bw,aw0,bw0,\n"
	  "                  #  ad,bdk,alpha,beta\n"
	  "  sampling hyperparameters:\n"
          "   -D cycles,start   #  sample alpha every this many cycles\n"
          "   -E cycles,start  #  sample beta every this many cycles\n"
	  "   -F var         #  fix var, var=a,b,a0,b0,aw,bw,aw0,bw0,ad,bdk,\n"
          "                  #           alpha, beta\n"
	  "   -g var,cnt     #  extra integer parameter for sampling var\n"
	  "   -G var,cycles,start #  sample var is same as -F\n"
	  "  control:\n"
          "   -c chkpnt      #  checkpoint all stats and pars every so many cycles\n"
          "   -C cycles      #  major Gibbs cycles\n"
	  "   -d dots        #  print a dot after this many docs\n"
	  "   -e             #  send error log to STDERR\n"
	  "   -f FMT         #  'ldac', 'witdit', 'docword', 'bag', 'lst' for data format\n"
          "   -I init,cycle,inc,free  #  controls constrains on topic changes\n"
	  "   -J init,cycle,prop,best #  do a merge step starting at init, every cycle\n"
	  "                  #     and ignore topics with proportion < prop (default=1/(T*100))\n"
	  "                  #     add the best merges at each stage (default=1)\n"
          "   -K topics      #  maximum number of topics\n"
	  "   -m             #  up the memory conservation by one\n"
	  "   -M maxtime     #  maximum training seconds (wall time), quit early if reached\n"
	  "   -N maxNwt,maxT #  maximum counts for Stirling number tables\n"
#ifdef H_THREADS
	  "   -q threads     #  set number of threads, default 1\n"
#endif
	  "   -r offset      #  restart using data from offset on (usually 0)\n"
          "                  #  load training statistics previously saved\n"
	  "   -r phi         #  means load '.phi' file\n"
	  "   -r theta       #  means load '.theta' file\n"
#ifdef EXPERIMENTAL
	  "   -r hdp         #  means load saved data from HDP,\n"
          "                  #    from 'mode-word-assignments.dat' file\n"
#endif
          "   -s seed        #  random number seed, default is a time value\n"
	  "   -v             #  up the verbosity by one\n"
#ifdef EXPERIMENTAL
	  "   -w size,incr,start  #  use a data window this big, \n"
	  "                  #     drift by incr after start cycles\n"
#endif
	  "   -W W           #  change max W\n"
	  "   -x             #  enable use of exclude topics for -Q\n"
	  "  testing and reports:\n"
	  "   -h HOLD,arg    #  use document completion in '-l' testing\n"
          "                  #  HOLD=dict, hold out words w with (w%%arg)==0\n"
          "                  #  HOLD=doc, hold out at place l with (l%%arg)==0\n"
          "                  #  HOLD=fract, hold out last fract words in doc\n"
	  "   -h all         #  no test set, done on training set\n"
          "   -l DIAG,cycles,start #  cycles for runtime calculations\n"
	  "                  #  DIAG is one of 'alpha','phi','prog','sparse',\n"
	  "                  #     'theta' or 'testprob' \n"
          "   -L DIAG,cycles,start #  cycles for diagnostic calculations\n"
#ifdef EXPERIMENTAL
	  "                  #  DIAG is one of 'class','like','lrs','query'\n"
#else
#ifdef QUERY
	  "                  #  DIAG is one of 'class','like','query'\n"
#else
	  "                  #  DIAG is one of 'class','like'\n"
#endif
#endif
          "   -o SC[,count]  #  SC=score type, 'cost', 'count', 'idf', 'Q', 'phi'\n"
          "                  #     optionally add number of words to print\n"
	  "   -O             #  report likelihood, not scaled perplexity\n"
	  "   -p             #  report coherency via PMI of topics\n"
          "                  #    use twice to set #topics to value set by -o\n"
	  "   -P secs        #  calc test perplexity every interval in secs\n"
#ifdef QUERY
	  "   -Q nres,file   #  do queries using query file, must use '-rphi'\n"
#endif
          "   -t traindocs   #  train documents, at start, default all-test\n"
          "   -T testdocs    #  test documents, at end, default 0\n"
	  "   -T TESTSTEM    #  or stem for a data file of same kind\n"
          "   -V             #  load vocab file to allow printing terms\n"
	  "   -X             #  do classifier results using data in STEM.class\n"
 	  );
}

static time_t wall_secs() {
  struct timeval tp;
  struct timezone tz;
  gettimeofday(&tp, &tz);
  return tp.tv_sec;
}

void *sampling_p(void *pargs)
{
  int i;
  float *p = fvec(ddN.T * 4);
  D_MiSi_t dD;  D_pargs_p *par =(D_pargs_p *) pargs;
  int procs = par->procs;
  clock_t t1 = clock();
  int start;
  int end;

  if ( par->window>0 ) {
     start = ddP.window_left - ddP.window_incr;
     end = ddP.window_right;
  } else {
    start = 0;
    if ( ddP.window )
      end = ddP.window;
    else
      end = ddN.DT;
  }
  
  if ( ddP.bdk!=NULL ) 
    misi_init(&ddM,&dD);

   /*
   *  sampling
   */
  par->thislp = 0;
  par->thisNd = 0;
  for (i=start+par->processid; i<end; i+=procs) {   
    int incremental;
    int usei = i % ddN.DT;
    if ( usei<0 ) usei += ddN.DT;
    if ( ddP.bdk!=NULL )  //WRAY ???
      misi_build(&dD,usei,0); 
    incremental = 0;
#ifdef EXPERIMENTAL
    if ( par->window ) {
      if ( i<ddP.window_left )
	incremental = -1;
      else if ( i>=ddP.window_right-ddP.window_incr )
	incremental = 1;
    }
#endif
    par->thislp += gibbs_lda(GibbsNone, par->Tmax, usei, ddD.NdT[usei], p, 
			     &dD, incremental, par->processid);
    par->thisNd += ddD.NdT[usei];
    if ( par->dots>0 && i>0 && (i%par->dots==0) ) yap_message(".");
    if ( ddP.bdk!=NULL )   //WRAY ???
      misi_unbuild(&dD,usei,0); 
  }
  free(p);
  if ( ddP.bdk!=NULL ) 
    misi_free(&dD);
  par->tot_time = (double)(clock() - t1) / CLOCKS_PER_SEC;
  return NULL;
}

void *testing_p(void *pargs)
{
  int i;
  float *p = fvec(ddN.T * 4);
  D_MiSi_t dD;
  D_pargs_p *par =(D_pargs_p *) pargs;
  int procs = par->procs;
  clock_t t1 = clock();
  enum GibbsType fix = par->fix;
  
  if ( ddP.bdk!=NULL ) 
    misi_init(&ddM,&dD);
  /*
   *  sampling
   */
  par->thislp = 0;
  par->thisNd = 0;
  for (i=ddN.DT+par->processid; i<ddN.D; i+=procs) {    
    int  thisw =  add_doc(i, fix);
    if ( thisw<=1 ) {
      remove_doc(i, fix);
      continue;
    }
    if ( ddP.bdk!=NULL ) misi_build(&dD,i,0);
    par->thislp += gibbs_lda(fix, par->Tmax, i, ddD.NdT[i], p, &dD, 
			     0, par->processid);
    par->thisNd += ddD.NdT[i];
    remove_doc(i, fix);
    if ( par->dots>0 && i>0 && (i%par->dots==0) ) 
      yap_message(".");
    if ( ddP.bdk!=NULL ) 
      misi_unbuild(&dD,i,0);
  }
  free(p);
  if ( ddP.bdk!=NULL ) 
    misi_free(&dD);
  par->tot_time = (double)(clock() - t1) / CLOCKS_PER_SEC;
  return NULL;
}

static float distH(float *p1, float *p2, int N) {
  int n;
  float tot = 0;
  for (n=0; n<N; n++) {
    float diff = sqrt(p1[n])-sqrt(p2[n]);
    tot += diff*diff;
  }
  tot /= N;
  return tot;
}

/*==========================================
 * main
 *========================================== */
int main(int argc, char* argv[])
{
  int c, iter, ITER=0;
  unsigned long seed=0;
  int Tmax = 0;
  enum dataType data = LdaC;
  enum dataType testdata = LdaC;
  int dots = 0;

  enum GibbsType fix_hold = GibbsNone;
  char *stem;
#ifdef QUERY
  int qparts=0, this_qpart=0;
#endif
  char *resstem;
  int displaycount = 10;
  int pmicount = 10;
  char *betafile = NULL;
  char *alphafile = NULL;
  enum PDPType PYalphain=H_None, PYbetain=H_None;
  double betacin = 0;
  double alphacin = 0;
  int noerrorlog = 0;
  double probepsilon = 0;
  int load_vocab = 0;
  int checkpoint = 0;
  int restart_offset=0;
  int restart = 0;
  int maxNwt = 10000;
  int maxT = 1000;
  double BDKval = 0;
  enum ScoreType score=ST_idf;
  
  int doexclude = 0;
  int cal_perp = 0;
  clock_t t1=0, t2=0, t3=0;
  time_t wall_start = 0;
  double tot_time = 0;
  double psample_time = 0;
  double t_interval = 0;
  long max_time = LONG_MAX;
  int num_perp = 1;
  enum ParType par;
  int procs = 1;          /*  number of threads */
  int dopmi = 0;
  int showlike = 0;
  int nosave = 0;
  int doclass = 0;
#ifdef EXPERIMENTAL
  int loadhdp = 0;
#endif
  int loadtheta = 0;
  int loadphi = 0;
  int nosample = 0;
  int maxW = 0;
#ifdef QUERY
  char *queryfile = NULL;
  int querycnt = 0;
#endif
  /*
   *  default values
   */
  ddN.T = 10;
  ITER = 100;
  ddN.TEST = 0;

  pctl_init();
  diag_alloc();

  while ( (c=getopt(argc, argv,"A:B:c:C:d:D:eE:f:F:g:G:h:iI:J:K:l:L:mM:N:o:OpP:q:Q:r:R:s:S:t:T:vVw:W:xX"))>=0 ) {
    switch ( c ) {
    case 'A':
      if ( !optarg )
	yap_quit("Need a valid 'A' argument\n");
      if ( strncmp(optarg,"hdp",3)==0 ) 
	ddP.PYalpha = H_HDP;
      else if ( strncmp(optarg,"hpdd",4)==0 ) 
	ddP.PYalpha = H_HPDD;
      else if ( strncmp(optarg,"ng",2)==0 ) 
	ddP.PYalpha = H_NG;
      else if ( strncmp(optarg,"pdp",3)==0 ) 
	ddP.PYalpha = H_PDP;
      else if ( strncmp(optarg,"dir",3)==0 ) 
	ddP.PYalpha = H_None;
      else if ( sscanf(optarg,"%lf",&alphacin)==1 )
	ddP.PYalpha = H_None;
      else
	yap_quit("Need a valid 'A' argument\n");
      /*  save copy so later reload of par file wont loose it */
      PYalphain = ddP.PYalpha;
      if ( ddP.PYalpha != H_HPDD && ddP.PYalpha != H_NG ) {
	/*  get file part */
	char *farg = strchr(optarg,',');
	if ( farg!=NULL )
	  alphafile = farg+1;
	if ( ddP.PYalpha==H_None && alphafile 
	     && strcmp(alphafile,"uniform")==0 )
	  alphafile = NULL;
      }
      break;
    case 'B':
      if ( !optarg )
	yap_quit("Need a valid 'B' argument\n");
      if ( strncmp(optarg,"hdp",3)==0 ) 
	ddP.PYbeta = H_HDP;
      else if ( strncmp(optarg,"hpdd",4)==0 ) 
	ddP.PYbeta = H_HPDD;
      else if ( strncmp(optarg,"pdp",3)==0 ) 
	ddP.PYbeta = H_PDP;
      else if ( strncmp(optarg,"dir",3)==0 ) 
	ddP.PYbeta = H_None;
      else if ( sscanf(optarg,"%lf",&betacin)==1 ) {
	ddP.PYbeta = H_None;
      } else
	yap_quit("Need a valid 'B' argument\n");
      /*  save copy so later reload of par file wont loose it */
      PYbetain = ddP.PYbeta;
      if ( ddP.PYbeta != H_HPDD ) {
	/*  get file part */
	char *farg = strchr(optarg,',');
	if ( farg!=NULL ) {
	  betafile = farg+1;
	  if ( ddP.PYbeta==H_None && betafile 
	       && strcmp(betafile,"uniform")==0 )
	    betafile = NULL;
	}
      }
      break;
    case 'c':
      if ( !optarg || sscanf(optarg,"%d",&checkpoint)!=1 )
        yap_quit("Need a valid 'c' argument\n");
      break;
    case 'C':
      if ( !optarg || sscanf(optarg,"%d",&ITER)!=1 )
	yap_quit("Need a valid 'C' argument\n");
      break;
    case 'd':
      if ( !optarg || sscanf(optarg,"%d",&dots)!=1 )
	yap_quit("Need a valid 'd' argument\n");
      break;
    case 'D':
      if ( !optarg || sscanf(optarg,"%d,%d",
			     &ddT[ParAlpha].cycles,&ddT[ParAlpha].start)<1 )
	yap_quit("Need a valid 'D' argument\n");
      break;
    case 'e':
      noerrorlog++;
      break;
   case 'E':
      if ( !optarg || sscanf(optarg,"%d,%d",
			     &ddT[ParBeta].cycles,&ddT[ParBeta].start)<1 )
        yap_quit("Need a valid 'E' argument\n");
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
      else if ( strcmp(optarg,"lst")==0 ) 
	data = SeqTxtBag;
       else
	yap_quit("Illegal data type for -f\n");
      break;
    case 'F':
      if ( strcmp(optarg,"all")==0 ) {  
	for (par=ParA; par<=ParBeta; par++) 
	  ddT[par].fix = 1;
	nosample = 1;
      } else if ( strcmp(optarg,"phi")==0 ) {
	loadphi++;
	/*
	 *  special case since we call these betatot and alphatot
	 */
      } else if ( strcmp(optarg,"beta")==0 ) {
	ddT[ParBeta].fix = 1;
      } else if ( strcmp(optarg,"alpha")==0 ) {
	ddT[ParAlpha].fix = 1;
      } else {
	par = findpar(optarg);
	if ( par==ParNone ) {
	  yap_quit("Illegal arg for -F\n");
	} else {
	  ddT[par].fix = 1;
	}
      }
      break;
    case 'g':
	{
	  char var[100];
	  int st=0;
	  if ( !optarg || sscanf(optarg,"%[^, ],%d", &var[0], &st)<1  )
            yap_quit("Need a valid 'g' argument\n");
          par = findpar(var);
          if ( par==ParBDK )
            ddP.kbatch = st;
          else
            yap_quit("Illegal var for -g\n");
        }
        break;      
    case 'G':
      {
	char var[100];
	int st=0, cy=0;
	if ( !optarg || sscanf(optarg,"%[^, ],%d,%d",
			       &var[0], &cy, &st)<2 || st<0 || cy<0 )
	  yap_quit("Need a valid 'G' argument\n");
	par = findpar(var);
	if ( par==ParNone )
	  yap_quit("Illegal var for -G\n");
	ddT[par].start = st;
	ddT[par].cycles = cy;
      }
      break;
    case 'h':
      {
	fix_hold = GibbsHold;
	if ( !optarg  )
	  yap_quit("Need a valid 'h' argument\n");
        if ( strcmp(optarg,"all")==0 ) {
          ddP.hold_all = 1;
	  ddP.mltburn = 0;
	  ddP.mltiter = 1;
        } else if ( strncmp(optarg,"dict,",5)==0 ) {
          if ( sscanf(&optarg[5],"%d",&ddP.hold_dict)<1 || ddP.hold_dict<2 )
            yap_quit("Need a valid 'hdict' argument\n");
        } else if ( strncmp(optarg,"fract,",6)==0 ) {
          if ( sscanf(&optarg[6],"%lf",&ddP.hold_fraction)<1 
               || ddP.hold_fraction<=0 || ddP.hold_fraction>=1 )
            yap_quit("Need a valid 'hfract' argument\n");
        } else if ( strncmp(optarg,"doc,",4)==0 ) {
          if ( sscanf(&optarg[4],"%d",&ddP.hold_every)<1 || ddP.hold_every<2 )
            yap_quit("Need a valid 'hdoc' argument\n");
        } else
          yap_quit("Need a valid 'h' argument\n");
      }
      break;
    case 'I':
      if ( !optarg || sscanf(optarg,"%d,%d,%d,%d",
                             &ddP.Tinit,&ddP.Tcycle,&ddP.Tinc,&ddP.Tfree)<1 )
        yap_quit("Need a valid 'I' argument\n");
      break;
    case 'J':
      if ( !optarg || sscanf(optarg,"%d,%d,%f,%d",
                             &ddP.mergeinit,&ddP.mergeiter,&ddP.mergemin, 
			     &ddP.mergebest)<1 
	   || ddP.mergeinit<1 || ddP.mergeiter<2 || ddP.mergemin>0.5 )
        yap_quit("Need a valid 'J' argument\n");
      break;
    case 'K':
      if ( !optarg || sscanf(optarg,"%d",&ddN.T)!=1 )
	yap_quit("Need a valid 'K' argument\n");
      break;
    case 'l':
      if ( !optarg )
	yap_quit("Need a valid 'l ' argument\n");
      if ( strncmp(optarg,"sparse,",7)==0 ) {
	if ( sscanf(&optarg[7],"%d,%d",&ddP.spiter, &ddP.spburn)<2 )
	  yap_quit("Need a valid 'l sparse,' argument\n");
      } else if ( strncmp(optarg,"testprob,",9)==0 ) {
	if ( sscanf(&optarg[9],"%d,%d,%lf",
		    &ddP.tprobiter, &ddP.tprobburn, &probepsilon)<2 )
	  yap_quit("Need a valid 'l testprob,' argument\n");
      } else if ( strncmp(optarg,"theta,",6)==0 ) {
	if ( sscanf(&optarg[6],"%d,%d,%lf",
		    &ddP.probiter, &ddP.probburn, &probepsilon)<2 )
	  yap_quit("Need a valid 'l theta,' argument\n");
      } else if ( strncmp(optarg,"prog,",5)==0 ) {
	if ( sscanf(&optarg[5],"%d,%d",&ddP.progiter, &ddP.progburn)<2 )
	  yap_quit("Need a valid 'l prog,' argument\n");
      } else if ( strncmp(optarg,"phi,",4)==0 ) {
	if ( sscanf(&optarg[4],"%d,%d",&ddP.phiiter, &ddP.phiburn)<2 )
	  yap_quit("Need a valid 'l word,' argument\n");
      } else if ( strncmp(optarg,"alpha,",6)==0 ) {
	if ( sscanf(&optarg[6],"%d,%d",&ddP.alphaiter, &ddP.alphaburn)<2 )
	  yap_quit("Need a valid 'l word,' argument\n");
      } else
	yap_quit("Need a valid DIAG code in 'l' argument\n");
      break;
    case 'L':
      if ( !optarg )
	yap_quit("Need a valid 'L ' argument\n");
      if ( strncmp(optarg,"class,",6)==0 ) {
	if ( sscanf(&optarg[6],"%d,%d",&ddP.prditer, &ddP.prdburn)<1 )
	  yap_quit("Need a valid 'L class,' argument\n");
      } else if ( strncmp(optarg,"like,",5)==0 ) {
	if ( sscanf(&optarg[5],"%d,%d",&ddP.mltiter, &ddP.mltburn)<1 )
	  yap_quit("Need a valid 'L like' argument\n");
      } else if ( strncmp(optarg,"lrs,",4)==0 ) {
	if ( sscanf(&optarg[4],"%d,%d",&ddP.lrsiter, &ddP.lrsburn)<1 )
	  yap_quit("Need a valid 'L lrs,' argument\n");
#ifdef QUERY
      } else if ( strncmp(optarg,"query,",6)==0 ) {
	if ( sscanf(&optarg[6],"%d",&ddP.queryiter)<1 )
	  yap_quit("Need a valid 'L query,' argument\n");
#endif
      } else
	yap_quit("Need a valid DIAG code in 'L' argument\n");
      break;
    case 'm':
      ddP.memory++;
      break;
    case 'M':
      if(!optarg || sscanf(optarg, "%ld", &max_time) != 1)
	yap_quit("Need a valid 'N' argument\n");
      break;
    case 'N':
      if ( !optarg || sscanf(optarg,"%d,%d", &maxNwt, &maxT)<1 )
	yap_quit("Need a valid 'N' argument\n");
      break;
     case 'o':
      {
        char *cp;
        if ( strncmp(optarg,"idf",3)==0 )
          score = ST_idf;
        else if ( strncmp(optarg,"phi",3)==0 )
          score = ST_phi;
        else if ( strncmp(optarg,"count",5)==0 )
          score = ST_count;
        else if ( strncmp(optarg,"Q",1)==0 )
          score = ST_Q;
        else if ( strncmp(optarg,"cost",4)==0 )
          score = ST_cost;
        else
          yap_quit("Need a valid parameter for 'o' argument\n");
        cp = strchr(optarg,',');
        if ( cp && sscanf(&cp[1],"%d",&displaycount)<1 )
          yap_quit("Need a second valid '-o' count argument\n");
      }
      break;
    case 'O':
      showlike++;
      break;
    case 'p':
      dopmi++;
      break;
    case 'P':
      cal_perp = 1;
      if(!optarg || sscanf(optarg, "%lf", &t_interval) != 1)
	yap_quit("Need a valid 'P' argument\n");
      break;
#ifdef H_THREADS
   case 'q':
      if(!optarg || sscanf(optarg, "%d", &procs) != 1)
	yap_quit("Need a valid 'q' argument\n");
      break;
#endif
#ifdef QUERY
    case 'Q':
      {
	queryfile = malloc(strlen(optarg));
	if( !optarg || (sscanf(optarg, "%d,%[^,],%d/%d", 
			       &querycnt, queryfile, &this_qpart, &qparts) < 2) )
	  yap_quit("Need a valid 'Q' argument\n"); 
	nosample = 1;
      }
      break;
#endif
    case 'r':
      restart++;
      if ( !optarg )
	yap_quit("Need a valid 'r' argument\n");
      if ( strcmp(optarg,"theta")==0 ) {
	loadtheta++;
      } else if ( strcmp(optarg,"phi")==0 ) {
	loadphi++;
#ifdef EXPERIMENTAL
      } else if ( strcmp(optarg,"hdp")==0 ) {
	loadhdp++;
	if ( ddP.PYbeta != H_None ) {
	  ddP.PYbeta = H_None;
	}
#endif
      } else if ( sscanf(optarg,"%d",&restart_offset)!=1 )
	yap_quit("Need a valid 'r' argument\n");
      break;
    case 'R':
      if ( !optarg || sscanf(optarg,"%d",&procs)!=1 )
	yap_quit("Need a valid 'R' argument\n");
      break;
    case 's':
      if ( !optarg || sscanf(optarg,"%lu",&seed)!=1 )
	yap_quit("Need a valid 's' argument\n");
      break;
    case 'S':
      {
	char var[100];
	double vin=0;
	if ( !optarg || sscanf(optarg,"%[^=, ]=%lf",
			       &var[0], &vin)<2  )
	  yap_quit("Need a valid 'S' argument\n");
	par = findpar(var);
	if ( par==ParNone )
	  yap_quit("Illegal var for -S\n");
	else if ( par==ParBDK ) {
	  BDKval = vin;
	} else
	  *(ddT[par].ptr) = vin;
      }   
      break;
    case 't':
      if ( !optarg || sscanf(optarg,"%d",&ddP.training)!=1 )
	yap_quit("Need a valid 't' argument\n");
      break;
    case 'T':
      if ( !optarg )
	yap_quit("Need a valid 'T' argument\n");
      {
	char *tname = data_name(optarg,data);
	FILE *fp = fopen(tname,"r");
	if ( fp==NULL ) {
	  free(tname);
	  tname = data_name(optarg,testdata);
	  fp = fopen(tname,"r");
        } else {
	  testdata = data;
        }
	free(tname);
	if ( fp!=NULL ) {
	  /*  its a valid test filename */
          ddP.teststem = optarg;
	  fclose(fp);
	} else if ( sscanf(optarg,"%d",&ddN.TEST)!=1 )
	  yap_quit("Need a valid 'T' argument\n");
      }
      break;
    case 'v':
      verbose++;
      break;
    case 'V':
      load_vocab++;
      break;
#ifdef EXPERIMENTAL
    case 'w':
      if ( !optarg || sscanf(optarg,"%d,%d,%d",&ddP.window, 
			     &ddP.window_incr, &ddP.window_cycle)<2 )
	yap_quit("Need a valid 'w' argument\n");
#endif
      break;
    case 'W':
      if ( !optarg || sscanf(optarg,"%d",&maxW)<1 )
	yap_quit("Need a valid 'W' argument\n");
      break;
    case 'x': 
      doexclude = 1;
      break;
    case 'X':
      doclass = 1;
      break;
    default:
      yap_quit("Unknown option '%c'\n", c);
    }
  }

  if (argc-optind != 2) {
    usage();
    exit(-1);
  }
  if ( optind>=argc ) {
    yap_quit("No arguments given\n");
  }
  stem = strdup(argv[optind++]);
  resstem = strdup(argv[optind++]);

  if ( dopmi && load_vocab==0 )
    load_vocab++;
  if ( dopmi>1 )
    pmicount = displaycount;

  if ( restart && maxW>0 ) 
    yap_quit("Cannot change dimension with option -W if restart is done\n");

  if ( noerrorlog==0 ) {
    char *wname = yap_makename(resstem, ".log");
    yap_file(wname);
    free(wname);
  }
  
  yap_commandline(argc, argv);
  yap_message("Version " HCA_VERSION ", ");
#ifdef H_THREADS
  yap_message("threads, ");
#endif
  yap_message(stype());

  if ( restart || loadphi ) {
    char buf[1000];
    char *fname = yap_makename(resstem,".par");
    FILE *fp = fopen(fname,"r");
    if ( !fp ) 
      yap_quit("Parameter file '%s' doesn't exist\n", fname);
    fclose(fp);
    free(fname);
    ddN.T = atoi(readpar(resstem,"T",buf,50));
    pctl_read(resstem, buf);
    /*   if command line had set to PDP, then restore */
    if ( PYalphain == H_PDP )
      ddP.PYalpha = H_PDP;
    if ( PYbetain == H_PDP )
      ddP.PYbeta = H_PDP;
     if ( ddP.training==0 ) {
      char *pv = readpar(resstem,"TRAIN",buf,50);
      if ( pv ) 
	 ddP.training = atoi(pv);
    } 
    if ( restart || maxW==0 ) {
      maxW = atoi(readpar(resstem,"W",buf,50));
    }
    if ( doexclude==0 ) {
      if ( ddP.n_excludetopic ) {
	ddP.n_excludetopic = 0;
	free(ddP.excludetopic);
	free(ddP.bits_et);
      }
    }
  }

  /*
   *   correct parameters after command line
   */
  if ( BDKval>0 ) {
    int t ;
    ddP.bdk = dvec(ddN.T);
    for (t=0; t<ddN.T; t++)
      ddP.bdk[t]  = BDKval;
  }
  pctl_fix(ITER, loadphi);
  pctl_samplereport();
#ifdef EXPERIMENTAL2
  Tmax = ddP.Tinit;
#else
  Tmax = ddN.T;
#endif
  
  assert(ddN.T>0);
  assert(ddN.TEST>=0);
  
  if ( !restart && ITER==0 )
    yap_quit("Zero iterations only allowed on restart\n");

  if ( loadphi && ddP.phiiter>0 )
    yap_quit("Options '-l phi,...' and '-r phi' incompatible\n");
  if ( loadtheta && ddP.probiter>0 )
    yap_quit("Options '-l theta,...' and '-r theta' incompatible\n");
  if ( loadphi && ddP.mergeiter>0 )
    yap_quit("Options '-J...' and '-r phi' incompatible\n");
  if ( loadtheta && ddP.mergeiter>0 )
    yap_quit("Options '-J...' and '-r theta' incompatible\n");
  if ( ddP.PYbeta && ddP.mergeiter>0 )
    yap_quit("Option '-J...' must have plain Dirichlet on Beta side\n");
  

  /*
   *   set random number generator
   */
   if ( seed ) {
    rng_seed(rngp,seed);
   } else {
    rng_time(rngp,&seed);
   }
   yap_message("Setting seed = %lu\n", seed);
  
  /*
   *  read data and get dimensions
   */
   {
     D_bag_t *dbp = data_read(stem, data);
     int training = pctl_training(dbp->D);
     if ( ddP.teststem ) {
       D_bag_t *dbpt = data_read(ddP.teststem, testdata);
       /* need to load a separate test set, strip to bare training */
       data_shrink(dbp, training);
       ddN.TEST = dbpt->D;
       data_append(dbp, dbpt);
       free(dbpt->w);  free(dbpt->d); free(dbpt);
     }
     if ( maxW>0 ) {
       if ( dbp->W <= maxW ) 
	 dbp->W = maxW;
       if ( dbp->W > maxW )
	 data_vocabshrink(dbp, maxW);
     }     
     /*
      *  transfer into system
      */
     ddN.D = dbp->D;
     ddN.W = dbp->W;
     ddN.N = dbp->N;
     ddN.NT = dbp->N;
     ddN.DT = training;
     ddD.w = dbp->w;
     ddD.d = dbp->d;
     free(dbp);
     if ( ddN.DT<ddN.D ) {
       /*  recompute NT */
       int i;
       for (i=0; i<ddN.N; i++)
	 if ( ddD.d[i]>=ddN.DT )
	   break;
       ddN.NT = i;
     }
   }

#ifdef EXPERIMENTAL
   if ( ddP.window && loadhdp )
     yap_quit("Option '-w' must not be used with loadhdp\n");
#endif
   if ( ddP.tprobiter>0 && ddN.TEST==0 )
     yap_quit("Option '-ltestprob,...' must have test data\n");
   
   if ( ddP.tprobiter>0 && fix_hold==GibbsHold  )
     /*
      *  the loop to do the test probs is done using usual
      *  test handling, which includes the hold-out handling
      */
     yap_quit("Option '-ltestprob,...' must not do hold out testing\n");
   
   /*
    *   finalise pctl stuff, needs to be called after dims set
    */
   if ( betacin>0 )
     ddP.betatot = betacin*ddN.W;
   if ( alphacin>0 )
     ddP.alphatot = alphacin*ddN.T;
   if ( ddP.mergeiter>0 ) {
     if ( ddP.mergemin==0 ) 
       ddP.mergemin = 0.01/ddN.T;
   }

   /*
    *   if ddP.betatot or ddP.alphatot ==0 they'll be set to defaults
    *   when PY=H_None
    */
   pctl_dims();
   if ( alphafile==NULL && (ddP.PYalpha==H_HDP||ddP.PYalpha==H_PDP) ) {
     if ( restart ) {
       /*
        *  use stored version of alpha
        */
       char *fname=yap_makename(resstem,".alpha");
       //   the NULL stops it from rewriting the file back
       pctl_fixalpha(fname, NULL);
       free(fname);
     } else {
       pctl_fixalpha("uniform", resstem);
     } 
   } else {
     pctl_fixalpha(alphafile, resstem);
   }
   if ( verbose && alphafile!=NULL && strcmp(alphafile,"uniform")!=0 ) {
     yap_message("Probability file for alpha prior is '%s'\n", alphafile);
   }
   if ( betafile==NULL && (ddP.PYbeta==H_HDP||ddP.PYbeta==H_PDP) ) {
     if ( restart ) {
       /*
        *  use stored version of beta
        */
       char *fname=yap_makename(resstem,".beta");
       //   the NULL stops it from rewriting the file back
       pctl_fixbeta(fname, NULL);
       free(fname);
     } else {
       pctl_fixbeta("uniform", resstem);
     } 
   } else {
     pctl_fixbeta(betafile, resstem);
   }
   if ( verbose && betafile!=NULL && strcmp(betafile,"uniform")!=0 ) {
     yap_message("Probability file for beta prior is '%s'\n", betafile);
   }

   if ( loadphi ) {
     phi_load(resstem);     
   } 
   if ( loadtheta ) {
     ddP.theta = fmat(ddN.D,ddN.T);
     prob_load(resstem,".theta",ddP.theta);
     prob_load(resstem,".testprob",&ddP.theta[ddN.DT]);
   }
   data_alloc();
   if ( ddP.phiiter>0 )
     phi_init(resstem);
   else 
     ddS.phi = NULL;
   if ( ddP.alphaiter>0 )
     alpha_init(resstem);
   else 
     ddS.alpha = NULL;
   if ( doclass )
     data_class(stem);
#ifdef QUERY
   if ( queryfile ) {
     if ( loadphi==0  ) 
       yap_quit("Querying with -Q needs phi (using -r phi) and no test\n");
     query_read(queryfile);
     if ( ddP.queryiter==0 )
       ddP.queryiter = 10;
   }
#endif
   hca_alloc();
   if ( ddP.bdk!=NULL ) 
     dmi_init(&ddM, ddS.z, ddD.w, ddD.NdTcum,
              ddN.T, ddN.N, ddN.W, ddN.D, ddN.DT,
              (fix_hold==GibbsHold)?pctl_hold:NULL);
   if ( load_vocab )
     data_vocab(stem);
   
   if ( probepsilon<=0 )
     probepsilon = 0.001/ddN.T;
   
   /*
    *    check if there are words to report sparsity on
    */
   if ( ddP.spiter>0 ) {
     char *fname = yap_makename(stem,".smap");
     FILE *fp = fopen(fname,"r");
     sparsemap_init(fp, procs);
     if ( fp ) 
       fclose(fp);
     free(fname);
   }
   if ( ddP.probiter>0 || ddP.tprobiter>0 ) {
     tprob_init();
   } 
   
   if ( ddP.phi==NULL && ddS.phi==NULL && score==ST_phi ) 
     yap_quit("Option '-o phi' needs a loaded phi matrix\n");

   /*
    *   setup the caches
    */
   cache_init(maxT, maxNwt);
   
   /*
    *  yap some details
    */
   data_report(ITER, seed);
   pctl_report();
   
  /*
   *  load/init topic assignments and prepare statistics
   */
#ifdef EXPERIMENTAL
  if ( loadhdp ) {
    hca_reset_stats(resstem, 0, 1, 0, ddN.DT);
    hca_load_hdp(resstem);
  } else 
#endif
  {
    if ( restart ) {
      hca_read_z(resstem, 0, ddN.DT);
      hca_rand_z(ddP.Tinit, ddN.DT, ddN.D);
    } else {
      hca_rand_z(ddP.Tinit, 0, ddN.D);
    }
    hca_reset_stats(resstem, restart, 0, 0, ddP.window?ddP.window:ddN.DT);
  }
  if ( ddP.PYalpha )
    yap_message("Initialised with %d classes\n", ddS.TDTnz);
  
#ifdef QUERY
  if ( queryfile ) {
    char *qname = yap_makename(queryfile,".out");
    gibbs_query(stem, querycnt,qname,dots,this_qpart,qparts);
    free(qname);
    restart_offset = ddN.DT;
    cal_perp = 0;
    checkpoint = 0;
    ddP.progiter = 1000;
    nosave= 1;
    goto FINISH;
  }
#endif

  if ( restart && ITER ) 
      yap_message("Initial log_2(perp)=%.4lf\n", likelihood() * (showlike?1:-M_LOG2E/ddN.NT));

  if ( ITER )
      yap_report("cycles: ");

  wall_start = wall_secs();
  
  for (iter=0; iter<ITER; iter++) {
    int pro;
    double thislp = 0;
    int   thisNd = 0;
    double testlp = 0;
    int   testNd = 0;  
#ifdef H_THREADS
    pthread_t thread[procs];
#endif
    D_pargs_p parg[procs];

    t1 = clock();
    
    /*
     *   set diagnostics to tell Gibbs which to store
     */
    /*  have to figure out what to do with threads */
    if ( ddG.n_words>0 && iter>ddP.spburn && (iter%ddP.spiter)==0 )
      ddG.docode = 1;
    else
      ddG.docode = 0;
    if ( ddP.probiter>0 && iter>ddP.probburn && (iter%ddP.probiter)==0 )
      ddG.doprob = 1;
    else
      ddG.doprob = 0;
    
#ifdef EXPERIMENTAL
    if ( ddP.window>0 && iter>=ddP.window_cycle ) {
      ddP.window_left += ddP.window_incr;
      ddP.window_right += ddP.window_incr;
      if ( ddP.window_left>= ddN.DT ) {
	ddP.window_left -= ddN.DT;
	ddP.window_right -= ddN.DT;
      }
    }
    //  yap_message("[ %d ... %d]\n", ddP.window_left,ddP.window_right);
#endif
    
    /*
     *  sampling
     */
#ifdef TRACE_WT
    for (i=0; i<ddN.NT; i++) {
      if ( ddD.w[i]==TR_W && Z_t(ddS.z[i])==TR_T )
	yap_message("Word=%d, topic=%d:  d=%d, z=%d, i=%d\n",
		    TR_W, TR_T, (int)ddD.d[i], (int)ddS.z[i], i);
    }
    yap_message("Word=%d, topic=%d:  start N=%d, T=%d\n",
		TR_W, TR_T, (int)ddS.Nwt[TR_W][TR_T],(int)ddS.Twt[TR_W][TR_T]);
#endif

    /*  a bit complex if no threads!  */
    for (pro = 0 ; pro < procs ; pro++){
      parg[pro].Tmax=Tmax;
      parg[pro].dots=dots;
      parg[pro].processid=pro;
      parg[pro].procs=procs;
      parg[pro].window = (iter>=ddP.window_cycle)?ddP.window:0;
#ifndef H_THREADS
      sampling_p(&parg[pro]);
#else
      if ( procs==1 )
          sampling_p(&parg[pro]);
      else if ( pthread_create(&thread[pro],NULL,sampling_p,(void*) &parg[pro]) != 0){
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

    // getting lp, Nd and clock
    for(pro = 0; pro < procs; pro++){
      thislp +=  parg[pro].thislp;
      thisNd +=  parg[pro].thisNd;
      tot_time += parg[pro].tot_time;
    }
#if 0 || defined(NONATOMIC)
    if ( procs>1 )
      hca_correct_twt();
#endif


    if ( ddG.docode ) {
      ddG.didcode++;
      ddG.docode = 0;
    }

    /*
     *  now run testing
     */
    if ( ddP.tprobiter>0 && iter>ddP.tprobburn && (iter%ddP.tprobiter)==0 ) {
      /*
       *   now run sampler on test docs
       */  
      ddG.doprob = 1;
      
      /*  a bit complex if no threads!  */
      for(pro = 0 ; pro < procs ; pro++){
        parg[pro].Tmax=Tmax;
        parg[pro].dots=dots;
        parg[pro].processid=pro;
        parg[pro].procs=procs;
#ifndef H_THREADS
        testing_p(&parg[pro]);
#else
        if ( procs==1 )
	    testing_p(&parg[pro]);
        else if( pthread_create(&thread[pro],NULL,testing_p,(void*) &parg[pro]) != 0){
          yap_message("thread failed %d\n",pro+1 );
        }
#endif
      }
#ifdef H_THREADS
      if ( procs>1 ) {
          //waiting for threads to finish
          for(pro = 0; pro < procs; pro++){
            pthread_join(thread[pro], NULL);
          }
      }
#endif 
      // getting lp, Nd and clock
      for(pro = 0; pro < procs; pro++){
        testlp +=  parg[pro].thislp;
        testNd +=  parg[pro].thisNd;
        tot_time += parg[pro].tot_time;
      }
      ddG.didtprob++;
      ddG.doprob = 0;
    }
    
    /*
     *   unset diagnostics
     */
    if ( ddG.doprob ) {
      ddG.didprob++;
      ddG.doprob = 0;
    }
    
    /*
     *   merge step
     */
    if ( ddP.mergeiter>0 && 
	 iter>ddP.mergeinit && (iter%ddP.mergeiter)==0 && iter-1<ITER ) {
      like_merge(ddP.mergemin, showlike?1:(-M_LOG2E/ddN.NT), ddP.mergebest);
    }

    /*
     *   sample hyperparameters
     */
    t3 = clock();
    if ( ddP.PYalpha ) {
      int t;
      for (t=0; t<ddN.T; t++)
        ddS.Tlife[t]++;
    }
    if ( nosample==0 ) 
      pctl_sample(iter,procs);
#ifdef EXPERIMENTAL2
    {
	int Tmax_before = Tmax;
    	Tmax = pctl_Tmax(Tmax, iter);
	if ( verbose>1 && Tmax_before < Tmax ) {
	  yap_message("T increased from %d to %d\n", Tmax_before, Tmax);
	}
    }
#endif
   
    /*
     *   do time calcs here to remove diagnostics+reporting
     */
    t2 = clock();
    tot_time += (double)(t2 - t1) / CLOCKS_PER_SEC;
    psample_time += (double)(t2 - t3) / CLOCKS_PER_SEC;
    if (cal_perp == 1) {
      if(tot_time >= (double)(t_interval * num_perp) || iter == 1){
	double logprob;
	char *teststr = fix_hold==GibbsHold?"Hold":"ML";
	t1 = clock();
	if ( ddP.window ) 
	  hca_reset_stats(resstem, 0, 0, 0, ddN.DT);
 	logprob = lp_test_ML(fix_hold, procs);
	if ( ddP.window ) 
	  hca_reset_stats(resstem, 0, 0, ddP.window_left,  ddP.window_right);
	t2 = clock();
	yap_message("\nTest: tot train/test time = %.4lf(s)/ %.4lf(s), "
                    "cycle = %d, log_2(test perp%s) = %.4f\n",
		    tot_time, (double)(t2 - t1) / CLOCKS_PER_SEC, iter, 
		    teststr, -M_LOG2E * logprob);
	num_perp++;
      }
    }
    /*
     *   progress reports
     */
    if ( ( iter>ddP.progburn && (iter%ddP.progiter)==0 ) || iter+1>=ITER ) {
      if ( ddP.window ) 
	hca_reset_stats(resstem, 0, 0, 0, ddN.DT);
      yap_message(" %d\nlog_2(perp)=%.4lf,%.4lf", iter, 
		  likelihood() * (showlike?1:-M_LOG2E/ddN.NT), thislp * (showlike?1:-M_LOG2E/thisNd));
      if ( ddP.tprobiter>0 ) 
	yap_message(", log_2(testperp)=%.4lf",
		    testlp * (showlike?1:-M_LOG2E/testNd));
      yap_message("\n");
      pctl_update(iter);
      if ( iter>0 && verbose>1 ) {
	if ( ddS.Ndt ) yap_probs();
	hca_displaytopics(stem, resstem, displaycount, score, 0,
			  (load_vocab>1)?1:0);
	if ( ddG.n_words>0 && ddG.didcode ) 
	  sparsemap_report(resstem,0.5,procs);
      }
      if ( ddP.window ) 
	hca_reset_stats(resstem, 0, 0, ddP.window_left,  ddP.window_right);
      if ( iter+1<ITER ) {
	yap_report("cycles: ");
      }
    } else {
      if ( verbose>1 )  
	yap_message("cycle %d\n", iter);
      else
	yap_message(" %d", iter);
    }
  
    if ( checkpoint>0 && iter>0 && iter%checkpoint==0 ) {
      data_checkpoint(resstem, stem, iter+1);
      yap_message(" checkpointed\n");
      hca_displaytopics(stem, resstem, displaycount, score, dopmi?pmicount:0,
			(load_vocab>1)?1:0);
      hca_report(resstem, stem, ITER, procs, fix_hold, showlike, nosave);
    }
    if ( ddP.phiiter>0 && iter>ddP.phiburn && (iter%ddP.phiiter)==0 )
      phi_update();
    if ( ddP.alphaiter>0 && iter>ddP.alphaburn && (iter%ddP.alphaiter)==0 )
      alpha_update();

    if ( wall_secs()-wall_start > max_time){ //72000 20hs, 86400 24 hs
      yap_message("\ntotal training time exceeds maximal training time  %ld s, quitting...\n", max_time);
      break;
    }

  } // over iter
  
  if ( ddP.window ) 
    hca_reset_stats(resstem, 0, 0, 0, ddN.DT);

  if ( ITER ) 
      yap_report("Finished after %d cycles on average of %lf+%lf(s) per cycle\n",
	     iter,  (tot_time-psample_time)/iter, psample_time/iter);
  
  if ( ( verbose==1 || ((iter+1)%5!=0 && verbose>1) ) ) {
    hca_displaytopics(stem, resstem, displaycount, score, dopmi?pmicount:0,
		      (load_vocab>1)?1:0);
    if ( ddG.n_words>0  && ddG.didcode) 
      sparsemap_report(resstem,0.5,procs);
  }
  
  if ( ddG.didtprob )
    tprob_report(resstem,probepsilon);
  if ( ddG.didprob )
    prob_report(resstem,probepsilon);

  if ( doclass ) {
    hca_displayclass(resstem);
  }
  
  if ( ddS.Ndt ) yap_probs();

  if ( ITER>0 && nosave==0 ) 
	data_checkpoint(resstem, stem, ITER);

  hca_report(resstem, stem, ITER, procs, fix_hold, showlike, nosave);
#ifdef QUERY
 FINISH:
#endif
  if ( ddP.phiiter>0 )
      phi_save(resstem);
  if ( ddP.alphaiter>0 )
      alpha_save(resstem);

  /*
   *  free
   */
  phi_free();
  alpha_free();
  pctl_free();
  cache_free();
  data_free();
  dmi_free(&ddM);
  hca_free();
  free(stem);
  free(resstem);
#ifdef QUERY
  if ( queryfile ) free(queryfile);
#endif
  rng_free(rngp);

  return 0;
}
