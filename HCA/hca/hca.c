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

void hca_displaytopics(char *resstem, int topword, enum ScoreType score);
void hca_displayclass(char *resstem);

//==================================================
// global variables
//==================================================

rngp_t rngp = NULL;
int verbose = 0;

#define QUERY

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
	"  H.Pitman-Yor sampler for topics"
	", H.Pitman-Yor sampler for words"
	"\n";
      else
	return 
	  "  H.Pitman-Yor sampler for topics"
	  ", Dirichlet sampler for words"
	  "\n";
  } else {
    if ( ddP.PYbeta )
      return 
	"  Dirichlet sampler for topics"
	", H.Pitman-Yor sampler for words"
	"\n";
    else
      return 
	"  Dirichlet sampler for topics"
	", Dirichlet sampler for words"
	"\n";
  }
}

static void usage() {
  fprintf(stderr,"Commandline:  OPTION+ STEM RESSTEM\n"
	  "  (reads STEM.dit and STEM.wit, results to RESSTEM.*)\n  ");
#ifdef H_THREADS
  fprintf(stderr," Threads,");
#endif
  fprintf(stderr,"%s", stype());
  fprintf(stderr,
          "  OPTION is choice of:\n"
	  "  setting hyperparameters:\n"
          "   -A value       #  use simple alpha prior and this value for it\n"
          "   -B value       #  use simple beta prior and this value for it\n"
          "   -A/B hdp/hpdd/pdp  #  for alpha or beta posn use this prior\n"
          "   -S var=value   #  initialise var=a,b,a0,b0,aw,bw,aw0,bw0,\n"
	  "                  #  ad,bdk\n"
	  "   -u bfile       #  beta proportions read from file\n"
	  "                  #     or use reserved words 'uniform' or 'file'\n"
	  "  sampling hyperparameters:\n"
          "   -D cycles,start   #  sample alpha every this many cycles\n"
          "   -E cycles,start  #  sample beta every this many cycles\n"
	  "   -F var        #  fix var, var=a,b,a0,b0,aw,bw,aw0,bw0,ad,bdk\n"
	  "   -g var,cnt     #  extra integer parameter for sampling var\n"
	  "   -G var,cycles,start #  sample var is same as -F\n"
	  "  control:\n"
          "   -c chkpnt      #  checkpoint all stats and pars every so many cycles\n"
          "   -C cycles      #  major Gibbs cycles\n"
	  "   -d dots        #  print a dot after this many docs\n"
	  "   -e             #  send error log to STDERR\n"
	  "   -f FMT         #  'ldac', 'witdit', 'docword', 'bag', 'lst' for data format\n"
          "   -I init,cycle,inc,free  #  controls constrains on topic changes\n"
          "   -K topics      #  maximum number of topics\n"
	  "   -m             #  up the memory conservation by one\n"
	  "   -M maxtime     #  maximum training time, quit early if reached\n"
	  "   -N maxT,maxNwt #  maximum counts for Stirling number tables\n"
#ifdef H_THREADS
	  "   -q threads     #  set number of threads, default 1\n"
#endif
	  "   -r offset      #  restart using data from offset on (usually 0)\n"
          "                  #  load training statistics previously saved\n"
	  "   -r alpha       #  keyword 'alpha' means load '.alpha' file\n"
	  "   -r phi         #  keyword 'phi' means load '.phi' file\n"
	  "   -r hdp         #  keyword 'hdp' means load saved data from HDP,\n"
          "                  #    from 'mode-word-assignments.dat' file\n"
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
          "   -l DIAG,cycles,start #  cycles for runtime calculations\n"
	  "                  #  DIAG is one of 'sp','theta','testprob','prog',\n"
	  "                  #  'phi', 'alpha'\n"
          "   -L DIAG,cycles,start #  cycles for diagnostic calculations\n"
#ifdef EXPERIMENTAL
	  "                  #  DIAG is one of 'lrs','class','like','query'\n"
#else
#ifdef QUERY
	  "                  #  DIAG is one of 'class','like','query'\n"
#else
	  "                  #  DIAG is one of 'class','like'\n"
#endif
#endif
          "   -o SC          #  SC=score type, 'count', 'idf', 'cost', 'Q', 'phi'\n"
	  "   -O             #  report likelihood, not scaled perplexity\n"
	  "   -p             #  report coherency via PMI of topics\n"
	  "   -P secs        #  calc test perplexity every interval in secs\n"
#ifdef QUERY
	  "   -Q nres,file   #  do queries using query file, must use '-rphi'\n"
#endif
          "   -t traindocs   #  train documents, at start, default all-test\n"
          "   -T testdocs    #  test documents, at end, default 0\n"
          "   -V             #  load vocab file to allow printing terms\n"
	  "   -X             #  do classifier results using data in STEM.class\n"
 	  );
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
			     &dD, incremental);
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
    par->thislp += gibbs_lda(fix, par->Tmax, i, ddD.NdT[i], p, &dD, 0);
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
  int qparts=0, this_qpart=0;
  char *resstem;
  char *betafile = NULL;
  int noerrorlog = 0;
  double probepsilon = 0;
  int displayed = 0;
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
  double tot_time = 0;
  double psample_time = 0;
  double t_interval = 0;
  double max_time = INFINITY;
  int num_perp = 1;
  enum ParType par;
  int procs = 1;          /*  number of threads */
  int dopmi = 0;
  int showlike = 0;
  int nosave = 0;
  int doclass = 0;
  int loadhdp = 0;
  int loadphi = 0;
  int loadalpha = 0;
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

  while ( (c=getopt(argc, argv,"A:B:c:C:d:D:eE:f:F:g:G:h:iI:K:l:L:mM:N:o:OpP:q:Q:r:R:s:S:t:T:u:vVw:W:xX"))>=0 ) {
    switch ( c ) {
    case 'A':
      if ( !optarg )
	yap_quit("Need a valid 'A' argument\n");
      if ( strcmp(optarg,"hdp")==0 ) 
	ddP.PYalpha = H_HDP;
      else if ( strcmp(optarg,"hpdd")==0 ) 
	ddP.PYalpha = H_HPDD;
      else if ( strcmp(optarg,"pdp")==0 ) 
	ddP.PYalpha = H_PDP;
      else if ( sscanf(optarg,"%lf",&ddP.alpha)==1 )
	ddP.PYalpha = H_None;
      else
	yap_quit("Need a valid 'A' argument\n");
      break;
    case 'B':
      if ( !optarg )
	yap_quit("Need a valid 'B' argument\n");
      if ( strcmp(optarg,"hdp")==0 ) 
	ddP.PYbeta = H_HDP;
      else if ( strcmp(optarg,"hpdd")==0 ) 
	ddP.PYbeta = H_HPDD;
      else if ( strcmp(optarg,"pdp")==0 ) 
	ddP.PYbeta = H_PDP;
      else if ( sscanf(optarg,"%lf",&ddP.beta)==1 )
	ddP.PYbeta = H_None;
      else
	yap_quit("Need a valid 'B' argument\n");
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
      } else if ( strcmp(optarg,"alpha")==0 ) {
	loadalpha++;
      } else if ( strcmp(optarg,"phi")==0 ) {
	loadphi++;
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
        if ( strncmp(optarg,"dict,",5)==0 ) {
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
    case 'K':
      if ( !optarg || sscanf(optarg,"%d",&ddN.T)!=1 )
	yap_quit("Need a valid 'K' argument\n");
      break;
    case 'l':
      if ( !optarg )
	yap_quit("Need a valid 'l ' argument\n");
      if ( strncmp(optarg,"sp,",3)==0 ) {
	if ( sscanf(&optarg[3],"%d,%d",&ddP.spiter, &ddP.spburn)<2 )
	  yap_quit("Need a valid 'l sp,' argument\n");
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
      if(!optarg || sscanf(optarg, "%lf", &max_time) != 1)
	yap_quit("Need a valid 's' argument\n");
      break;
    case 'N':
      if ( !optarg || sscanf(optarg,"%d,%d", &maxT, &maxNwt)<1 )
	yap_quit("Need a valid 'N' argument\n");
      break;
     case 'o':
      {
        if ( strcmp(optarg,"idf")==0 )
          score = ST_idf;
        else if ( strcmp(optarg,"phi")==0 )
          score = ST_phi;
        else if ( strcmp(optarg,"count")==0 )
          score = ST_count;
        else if ( strcmp(optarg,"Q")==0 )
          score = ST_Q;
        else if ( strcmp(optarg,"cost")==0 )
          score = ST_cost;
        else
          yap_quit("Need a valid parameter for 'o' argument\n");
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
      if ( strcmp(optarg,"alpha")==0 ) {
	loadalpha++;
      } else if ( strcmp(optarg,"phi")==0 ) {
	loadphi++;
      } else if ( strcmp(optarg,"hdp")==0 ) {
	loadhdp++;
	if ( ddP.PYbeta != H_None ) {
	  ddP.PYbeta = H_None;
	  ddP.beta = 1;
	}
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
	  tname = data_name(optarg,testdata);
	  fp = fopen(tname,"r");
        } else {
	  testdata = data;
        }
	if ( fp!=NULL ) {
	  /*  its a valid test filename */
          ddP.teststem = optarg;
	  fclose(fp);
	} else if ( sscanf(optarg,"%d",&ddN.TEST)!=1 )
	  yap_quit("Need a valid 'T' argument\n");
      }
      break;
    case 'u':
      ddP.betac = 0;
      betafile = optarg;
      break;
    case 'v':
      verbose++;
      break;
    case 'V':
      load_vocab = 1;
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

  if ( dopmi )
    load_vocab = 1;

  if ( loadalpha && loadphi==0 ) 
    yap_quit("If using -ralpha/-Falpha, should use -rphi/-Fphi\n");

  if ( noerrorlog==0 ) {
    char *wname = yap_makename(resstem, ".log");
    yap_file(wname);
    free(wname);
  }
  
  yap_commandline(argc, argv);
#ifdef H_THREADS
  yap_message(" Threads,");
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
    if ( ddP.training==0 ) {
      char *pv = readpar(resstem,"TRAIN",buf,50);
      if ( pv ) 
	 ddP.training = atoi(pv);
    } 
   if ( maxW==0 )
      maxW = atoi(readpar(resstem,"W",buf,50));
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
  if ( ddP.PYbeta==H_PDP && betafile==NULL )
    betafile = "uniform";
  if ( restart && betafile==NULL && (ddP.PYbeta==H_HDP||ddP.PYbeta==H_PDP)) 
    pctl_fix("restart", ITER);
  else
    pctl_fix(betafile, ITER);
  pctl_samplereport();
  Tmax = ddP.Tinit;
  
  assert(ddN.T>0);
  assert(ddN.TEST>=0);
  assert(restart || ITER>0);
  if ( loadphi && ddP.phiiter>0 )
    yap_quit("Options '-l phi,...' and '-r phi' incompatible\n");
  if ( loadalpha && ddP.alphaiter>0 )
    yap_quit("Options '-l alpha,...' and '-r alpha' incompatible\n");

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
   
   if ( ddP.window && loadhdp )
     yap_quit("Option '-w' must not be used with loadhdp\n");
   
   if ( ddP.tprobiter>0 && ddN.TEST==0 )
     yap_quit("Option '-ltestprob,...' must have test data\n");
   
   if ( ddP.tprobiter>0 && fix_hold==GibbsHold  )
     /*
      *  the loop to do the test probs is done using usual
      *  test handling, which includes the hold-out handling
      */
     yap_quit("Option '-ltestprob,...' must not do hold out testing\n");
	
  /*
   *   all data structures
   */
  pctl_dims();
  if ( loadphi ) {
    if ( score!=ST_phi ) {
      yap_message("Setting scoring to be 'phi' due to '-r phi' option\n");
      score = ST_phi;
    }
    phi_load(resstem);
    if ( loadalpha ) {
      alpha_load(resstem);
    }
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
  if ( load_vocab ) {
    data_vocab(stem);
  }

  if ( probepsilon<=0 )
    probepsilon = 0.001/ddN.T;

  /*
   *    check if there are words to report sparsity on
   */
  if ( ddP.spiter>0 ) {
    char *fname = yap_makename(stem,".smap");
    FILE *fp = fopen(fname,"r");
    if ( fp ) {
      sparsemap_init(fp);
      fclose(fp);
    } 
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
  if ( restart && betafile==NULL && (ddP.PYbeta==H_HDP||ddP.PYbeta==H_PDP)) {
    char *fname=yap_makename(resstem,".beta");
    //   the NULL stops it from rewriting the file back
    fixbeta(fname, NULL);
    free(fname);
  } else {
    fixbeta(betafile, resstem);
  }
  cache_init(maxT, maxNwt);
  
  /*
   *  yap some details
   */
  data_report(ITER, seed);
  pctl_report();
 
  /*
   *  load/init topic assignments and prepare statistics
   */
  if ( loadhdp ) {
    hca_reset_stats(resstem, 0, 1, 0, ddN.DT);
    hca_load_hdp(resstem);
  } else {
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
    data_df(stem);
    gibbs_query(querycnt,qname,dots,this_qpart,qparts);
    free(qname);
    restart_offset = ddN.DT;
    cal_perp = 0;
    checkpoint = 0;
    ddP.progiter = 1000;
    nosave= 1;
    // goto FINISH;
  }
#endif

  if ( restart && ITER ) 
      yap_message("Initial log_2(perp)=%.4lf\n", likelihood() * (showlike?1:-M_LOG2E/ddN.NT));

  if ( ITER )
      yap_report("cycles: ");
  
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
#ifndef H_THREADS
    /*  have to figure out what to do with threads */
    if ( ddG.n_words>0 && iter>ddP.spburn && (iter%ddP.spiter)==0 )
      ddG.docode = 1;
    else
      ddG.docode = 0;
    if ( ddP.probiter>0 && iter>ddP.probburn && (iter%ddP.probiter)==0 )
      ddG.doprob = 1;
    else
      ddG.doprob = 0;
#endif
    
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
      //calling sampling for processes
      if( pthread_create(&thread[pro],NULL,sampling_p,(void*) &parg[pro]) != 0){
        yap_message("thread failed %d\n",pro+1 );
      }
#endif
    }
#ifdef H_THREADS
    //waiting for threads to finish
    for (pro = 0; pro < procs; pro++){
      pthread_join(thread[pro], NULL);
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


#ifndef H_THREADS
    if ( ddG.docode ) {
      ddG.didcode++;
      ddG.docode = 0;
    }
    /*  have to figure out what to do with threads */
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
#ifdef NO_THREADS
        testing_p(&parg[pro]);
#else
        //calling sampling for processes
        if( pthread_create(&thread[pro],NULL,testing_p,(void*) &parg[pro]) != 0){
          yap_message("thread failed %d\n",pro+1 );
        }
#endif
      }
#ifndef NO_THREADS
      //waiting for threads to finish
      for(pro = 0; pro < procs; pro++){
        pthread_join(thread[pro], NULL);
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
#endif
    
    /*
     *   unset diagnostics
     */
    if ( ddG.doprob ) {
      ddG.didprob++;
      ddG.doprob = 0;
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
    {
	int Tmax_before = Tmax;
    	Tmax = pctl_Tmax(Tmax, iter);
	if ( verbose>1 && Tmax_before < Tmax )
		yap_message("T increased from %d to %d\n", Tmax_before, Tmax);
    }
   
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
      if ( verbose && iter%10==0 )
	yap_probs();
      if ( iter>0 && verbose>1 ) {
	if ( ddN.tokens ) {
            hca_displaytopics(resstem,20,score);
	    displayed++;
	}
	if ( ddG.n_words>0 && ddG.didcode ) 
	  sparsemap_report(resstem,0.5);
      }
      if ( ddP.window ) 
	hca_reset_stats(resstem, 0, 0, ddP.window_left,  ddP.window_right);
      if ( iter+1<ITER ) {
	yap_report("cycles: ");
      }
    } else {
      yap_message(" %d", iter);
      if ( verbose>1 )  yap_message("\n");
    }
  
    if ( checkpoint>0 && iter>0 && iter%checkpoint==0 ) {
      data_checkpoint(resstem, stem, iter+1);
      yap_message(" checkpointed\n");
      hca_report(resstem, stem, ITER, procs, fix_hold, 
		 (dopmi&&displayed>0)?1:0, showlike, nosave);
    }
    if ( ddP.phiiter>0 && iter>ddP.phiburn && (iter%ddP.phiiter)==0 )
      phi_update();
    if ( ddP.alphaiter>0 && iter>ddP.alphaburn && (iter%ddP.alphaiter)==0 )
      alpha_update();

    if (tot_time > max_time){ //72000 20hs, 86400 24 hs
      yap_message("\ntotal training time exceeds maximal training time  %lf s, quitting...\n", max_time);
      break;
    }

  } // over iter
  
  if ( ddP.window ) 
    hca_reset_stats(resstem, 0, 0, 0, ddN.DT);

  if ( ITER ) 
      yap_report("Finished after %d cycles on average of %lf+%lf(s) per cycle\n",
	     iter,  (tot_time-psample_time)/iter, psample_time/iter);
  
  if ( ( verbose==1 || ((iter+1)%5!=0 && verbose>1) ) ) {
    if ( ddN.tokens ) {
       hca_displaytopics(resstem,20,score);
       displayed++;
    }
    if ( ddG.n_words>0  && ddG.didcode) 
      sparsemap_report(resstem,0.5);
  }

  if ( ddG.didtprob )
    tprob_report(resstem,probepsilon);
  if ( ddG.didprob )
    prob_report(resstem,probepsilon);

  if ( doclass ) {
    hca_displayclass(resstem);
  }
  
  yap_probs();

  if ( ITER>0 && nosave==0 ) 
	data_checkpoint(resstem, stem, ITER);
  /*
   *   must have run hca_displaytopics() to create top words for PMI
   */
  hca_report(resstem, stem, ITER, procs, fix_hold, (dopmi&&displayed>0)?1:0, 
             showlike, nosave);

 FINISH:
  /*
   *  free
   */
  if ( ddP.phiiter>0 )
      phi_save(resstem);
  if ( ddP.alphaiter>0 )
      alpha_save(resstem);
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
