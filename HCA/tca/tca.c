/**
 * Main driver
 * Copyright (C) 2011-2014 Wray Buntine
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
#include "srng.h"
#include "stable.h"
#include "lgamma.h"
#include "tca.h"
#include "pctl.h"
#include "data.h"
#include "misi.h"
#include "stats.h"
#include "sample.h"
#include "probs.h"

/*
 *  not always threading, but code uses this anyway
 */
#include "pargs.h"

#ifdef H_THREADS
#include <pthread.h>
#endif
#include "atomic.h"

void tca_displaytopics(char *resstem, int topword, enum ScoreType score);
void checkm_evt(int w, int val);

//==================================================
// global variables
//==================================================

rngp_t rngp = NULL;
int verbose = 0;

/*
 *    Dimensions
 */
D_dims_t ddN;
/*
 *  hyperparameters
 */
D_pars_t ddP;
D_pctl_t ddT[ParBB+1];

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
 *   bursty data stuctures
 */
D_DMi_t ddM;

static void usage() {
  fprintf(stderr,"Commandline:  OPTION+ STEM RESSTEM\n"
	  "  (reads STEM.dit and STEM.wit, results to RESSTEM.*)\n  ");
  fprintf(stderr,
          "  OPTION is choice of:\n"
	  "  setting hyperparameters:\n"
          "   -S var=value   #  initialise var=aX,bXN for X=M,P and N=0,1\n"
	  "                  #          or var=aX,bX for X=T,B\n"
	  "  sampling hyperparameters:\n"
	  "   -F var        #  fix var, var as in -S\n"
	  "   -G var,cycles,start #  sample var is same as -S\n"
	  "  control:\n"
          "   -c chkpnt      #  checkpoint z and pars every so many cycles\n"
          "   -C cycles      #  major Gibbs cycles\n"
	  "   -d dots        #  print a dot after this many docs\n"
	  "   -e             #  send error log to STDERR\n"
	  "   -f FMT         #  'ldac', 'witdit', 'docword', 'bag' for data format\n"
          "   -K topics      #  maximum number of topics\n"
	  "   -N maxN,maxM   #  maximum counts for Stirling number tables\n"
	  "                  #     maxM is max #tables for all\n"
	  "                  #     maxN is max count for a_mu and a_phi\n"
	  "   -q threads     #  set number of threads, default 1\n"
	  "   -r             #  restart using data saved\n"
	  "   -R             #  restart from hca\n"
          "   -s seed        #  random number seed, default is a time value\n"
	  "   -v             #  up the verbosity by one\n"
	  "   -W W           #  make max W larger\n"
	  "  testing and reports:\n"
	  "   -h HOLD,arg    #  use document completion in '-l' testing\n"
          "                  #  HOLD=dict, hold out words w with (w%%arg)==0\n"
          "                  #  HOLD=doc, hold out at place l with (l%%arg)==0\n"
          "                  #  HOLD=fract, hold out last fract words in doc\n"
          "   -l DIAG,cycles,burn #  cycles for runtime calculations\n"
	  "                  #  DIAG is one of 'prog'\n"
          "   -L DIAG,cycles,burn #  cycles for diagnostic calculations\n"
	  "                  #  DIAG is one of 'like'\n"
          "   -o SC          #  SC=score type, 'count', 'idf', 'cost', 'Q'\n"
          "   -t traindocs   #  train documents, at start, default all-test\n"
          "   -T testdocs    #  test documents, at end, default 0\n"
          "   -V             #  load vocab file to allow printing terms\n"
 	  );
}

void *sampling_p(void *pargs)
{
  int i;
  float *p = fvec(ddN.T * 4);
  D_MiSi_t dD;
  D_pargs_p *par =(D_pargs_p *) pargs;
  int procs = par->procs;
  clock_t t1 = clock();
  
  if ( PCTL_BURSTY() )
    misi_init(&ddM,&dD);
  /*
   *  sampling
   */
  par->thislp = 0;
  par->thisNd = 0;
  for (i=par->processid; i<ddN.DT; i+=procs) {    
    if ( PCTL_BURSTY() )
      misi_build(&dD,i,0);
    par->thislp += gibbs_lda(GibbsNone, i, ddD.N_dT[i], p, &dD);
    par->thisNd += ddD.N_dT[i];
    if ( par->dots>0 && i>0 && (i%par->dots==0) ) 
      yap_message(".");
    if ( PCTL_BURSTY() )
      misi_unbuild(&dD,i,0);
  }
  free(p);
  if ( PCTL_BURSTY() )
    misi_free(&dD);
  par->tot_time = (double)(clock() - t1) / CLOCKS_PER_SEC;
  return NULL;
}

/*==========================================
 * main
 *========================================== */
int main(int argc, char* argv[])
{
  int c, iter, ITER=0, seed=0;
  enum dataType data = LdaC;
  enum dataType testdata = LdaC;
  int dots = 0;

  enum GibbsType fix_hold = GibbsNone;
  char *stem;
  char *resstem;
  int noerrorlog = 0;
  int displayed = 0;
  int load_vocab = 0;
  int checkpoint = 0;
  int restart = 0;
  int restart_hca = 0;
  int procs = 1;
  int maxW = 0;
  enum ScoreType score=ST_idf;
  
  double BM0val=0, BM1val =0, BP0val=0, BP1val=0;
  
  clock_t t1=0, t2=0, t3=0;
  double tot_time = 0;
  double psample_time = 0;
  enum ParType par;
  /*
   *  default values
   */
  ddN.T = 10;
  ITER = 100;
  ddN.TEST = 0;

  pctl_init();

  while ( (c=getopt(argc, argv,"c:C:d:ef:F:G:h:K:l:L:N:o:P:q:vrRs:S:t:T:vVW:"))>=0 ) {
    switch ( c ) {
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
    case 'e':
      noerrorlog++;
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
	for (par=ParAM; par<=ParBB; par++) 
	  ddT[par].fix = 1;
      } else {
	par = findpar(optarg);
	if ( par==ParNone )
	  yap_quit("Illegal arg for -F\n");
	ddT[par].fix = 1;
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
	if ( par==ParNone || par==ParB0P || par==ParB0M )
	  yap_quit("Illegal var for -G\n");
        ddT[par].fix = 0;
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
   case 'K':
      if ( !optarg || sscanf(optarg,"%d",&ddN.T)!=1 )
	yap_quit("Need a valid 'K' argument\n");
      break;
    case 'l':
      if ( !optarg )
	yap_quit("Need a valid 'l ' argument\n");
      if ( strncmp(optarg,"prog,",5)==0 ) {
	if ( sscanf(&optarg[5],"%d,%d",&ddP.progiter, &ddP.progburn)<2 )
	  yap_quit("Need a valid 'l prog,' argument\n");
      } else
	yap_quit("Need a valid DIAG code in 'l' argument\n");
      break;
    case 'L':
      if ( !optarg )
	yap_quit("Need a valid 'L ' argument\n");
      if ( strncmp(optarg,"like,",5)==0 ) {
	if ( sscanf(&optarg[5],"%d,%d",&ddP.mltiter, &ddP.mltburn)<1 )
	  yap_quit("Need a valid 'L like' argument\n");
      } else
	yap_quit("Need a valid DIAG code in 'L' argument\n");
      break;
    case 'N':
      if ( !optarg || sscanf(optarg,"%d,%d", &ddP.maxN, &ddP.maxM)<1 )
	yap_quit("Need a valid 'N' argument\n");
      break;
    case 'o':
      {
        if ( strcmp(optarg,"idf")==0 )
          score = ST_idf;
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
   case 'q':
      if(!optarg || sscanf(optarg, "%d", &procs) != 1)
	yap_quit("Need a valid 'q' argument\n");
      break;
    case 'r':
      restart++;
      break;
    case 'R':
      restart_hca++;
      break;
    case 's':
      if ( !optarg || sscanf(optarg,"%d",&seed)!=1 )
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
	else if ( par==ParBM0 ) 
	  BM0val = vin;
	else if ( par==ParBM1 ) 
	  BM1val = vin;
	else if ( par==ParBP0 ) 
	  BP0val = vin;
	else if ( par==ParBP1 ) 
	  BP1val = vin;
	else
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
    case 'v':
      verbose++;
      break;
    case 'V':
      load_vocab = 1;
      break;
    case 'W':
      if ( !optarg || sscanf(optarg,"%d",&maxW)<1 )
	yap_quit("Need a valid 'W' argument\n");
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

  if ( noerrorlog==0 ) {
    char *wname = yap_makename(resstem, ".log");
    yap_file(wname);
    free(wname);
  }
  
  yap_commandline(argc, argv);
#ifdef H_THREADS
  yap_message(" Threads,");
#endif

  if ( restart || restart_hca ) {
    char buf[1000];
    char *fname = yap_makename(resstem,".par");
    FILE *fp = fopen(fname,"r");
    if ( !fp ) 
      yap_quit("Parameter file '%s' doesn't exist\n", fname);
    fclose(fp);
    free(fname);
    ddN.T = atoi(readpar(resstem,"T",buf,50));
    if ( restart ) {
      ddN.E = atoi(readpar(resstem,"E",buf,50));
      pctl_read(resstem, buf);
    }
    if ( ddP.training==0 ) {
      char *pv = readpar(resstem,"TRAIN",buf,50);
      if ( pv ) 
	ddP.training = atoi(pv);
    } 
    if ( ddN.TEST==0 )
      ddN.TEST = atoi(readpar(resstem,"TEST",buf,50));
  } 

  assert(ddN.T>0);
  assert(ddN.TEST>=0);
  assert(restart || restart_hca || ITER>0);
	
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
    /*
     *  transfer into system
     */
    ddN.D = dbp->D;
    ddN.W = dbp->W;
    if ( ddN.W< maxW ) 
      ddN.W = maxW;
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

  data_read_epoch(stem);

  /*
   *   correct parameters after command line
   */
  pctl_fix(ITER);
  if ( BM0val>0 ) {
    ddP.b_mu[0] = BM0val;
  }
  if ( BM1val>0 ) {
    int i;
    for (i=1; i<ddN.E; i++)
      ddP.b_mu[i] = BM1val;
  }
  if ( BP0val>0 ) {
    int i;
    for (i=0; i<ddN.T; i++)
      ddP.b_phi[0][i] = BP0val;
  }
  if ( BP1val>0 ) {
    int i;
    if ( ddN.E==1 )
      yap_quit("b_phi[1] invalid when epochs==1\n");
    for (i=0; i<ddN.T; i++)
      ddP.b_phi[1][i] = BP1val;
  }
  pctl_samplereport();

  /*
   *   all data structures
   */
  data_alloc();
  if ( PCTL_BURSTY() ) 
    dmi_init(&ddM, ddS.z, ddD.w, ddD.N_dTcum,
             ddN.T, ddN.N, ddN.W, ddN.D, ddN.DT,
	     (fix_hold==GibbsHold)?pctl_hold:NULL);
  tca_alloc();
  if ( load_vocab ) {
    data_vocab(stem);
  }

  cache_init();
  
  /*
   *  yap some details
   */
  data_report(ITER, seed);
  pctl_report();
 
  /*
   *  load/init topic assignments and prepare statistics
   */
  if ( restart || restart_hca) {
    tca_read_z(resstem, 0, ddN.DT);
    tca_rand_z(ddN.T, ddN.DT, ddN.D);
  } else {
    tca_rand_z(ddN.T, 0, ddN.D);
  }
  tca_reset_stats(resstem, restart);

  if ( (restart || restart_hca ) && ITER ) 
      yap_message("Initial log_2(perp)=%lf\n", -M_LOG2E * likelihood()/ddN.NT);

  if ( ITER )
      yap_report("cycles: ");
  
  for (iter=0; iter<ITER; iter++) {
    int  pro;
    double thislp = 0;
    int   thisNd = 0;
#ifdef H_THREADS
    pthread_t thread[procs];
#endif
    D_pargs_p parg[procs];

    t1 = clock();
    
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
		TR_W, TR_T, (int)ddS.m_evt[TR_W][TR_T],(int)ddS.s_evt[TR_W][TR_T]);
#endif
#ifdef IND_STATS
    ddP.doc_ind_stats = u32tri(ddN.T,ddN.E,ddN.E);
    ddP.word_ind_stats = u32tri(ddN.T,ddN.E,ddN.E);
#endif

   /*  a bit complex if no threads!  */
    for (pro = 0 ; pro < procs ; pro++){
      parg[pro].dots=dots;
      parg[pro].processid=pro;
      parg[pro].procs=procs;
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
#ifdef H_THREADS
    if ( procs>1 )
      tca_reset_stats(NULL,1);
#endif

#ifdef IND_STATS
    {
      char *fname = yap_makename(resstem,".istats");
      FILE *ifp = fopen(fname,"a");
      int e1, e2, kk;
      fprintf(ifp,"Iteration %d\n", iter);
      for (kk=0; kk<ddN.T; kk++) {
	fprintf(ifp," Topic %d\n", kk);
	for (e1=0; e1<ddN.E; e1++) {
	  fprintf(ifp,"  Epoch %d\n     ", e1);
	  for (e2=0; e2<ddN.E; e2++)
	    fprintf(ifp," %u", (unsigned)ddP.doc_ind_stats[kk][e1][e2]);
	  fprintf(ifp,"\n     ");
	  for (e2=0; e2<ddN.E; e2++)
	    fprintf(ifp," %u", (unsigned)ddP.word_ind_stats[kk][e1][e2]);
	  fprintf(ifp,"\n");
	}
      }
      fclose(ifp);
      free(ddP.doc_ind_stats[0][0]); free(ddP.doc_ind_stats[0]); 
      free(ddP.doc_ind_stats); 
      free(ddP.word_ind_stats[0][0]); free(ddP.word_ind_stats[0]); 
      free(ddP.word_ind_stats);
      free(fname);
    }
#endif
    
    /*
     *   sample hyperparameters
     */
    t3 = clock();
    pctl_sample(iter);
   
    /*
     *   do time calcs here to remove diagnostics+reporting
     */
    t2 = clock();
    tot_time += (double)(t2 - t1) / CLOCKS_PER_SEC;
    psample_time += (double)(t2 - t3) / CLOCKS_PER_SEC;
    /*
     *   progress reports
     */
    if ( ( iter>ddP.progburn && (iter%ddP.progiter)==0 ) || iter+1>=ITER ) {
      yap_message("\n%d, log_2(perp)=%lf,%lf", iter, 
		  -M_LOG2E * likelihood()/ddN.NT, -M_LOG2E * thislp/thisNd);
      pctl_update(iter);
      if ( verbose && iter%10==0 )
	yap_probs();
      if ( iter>0 && verbose>1 ) {
	if ( ddN.tokens )
            tca_displaytopics(resstem,20,score);
	displayed++;
      }
      if ( iter+1<ITER ) {
	// yap_message("\n");
	yap_report("cycles: ");
      }
    } else {
      yap_message(" %d", iter);
      if ( verbose )  yap_message("\n");
    }
  
    if ( checkpoint>0 && iter>0 && iter%checkpoint==0 ) {
      data_checkpoint(resstem, stem, iter+1);
      yap_message(" checkpointed\n");
      tca_report(resstem, stem, ITER, procs, fix_hold);
    }

  } // over iter
  
  if ( ITER ) 
      yap_report("Finished after %d cycles on average of %lf+%lf(s) per cycle\n",
	     iter,  (tot_time-psample_time)/iter, psample_time/iter);
  
  if ( ( verbose==1 || ((iter+1)%5!=0 && verbose>1) ) ) {
    displayed++;
    if ( ddN.tokens )
       tca_displaytopics(resstem,20,score);
  }

  yap_probs();

  if ( ITER>0 ) 
	data_checkpoint(resstem, stem, ITER);
 
  tca_report(resstem, stem, ITER, procs, fix_hold);

  /*
   *  free
   */
  cache_free();
  pctl_free();
  data_free();
  dmi_free(&ddM);
  tca_free();
  free(stem);
  free(resstem);
  rng_free(rngp);

  return 0;
}
