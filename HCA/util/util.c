/*
 * Utilities
 * Copyright (C) 2010 Dave Newman
 * Copyright (C) 2011-2014 Wray Buntine
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Dave Newman (dave.newman@nicta.com.au)
 *         Wray Buntine (wray.buntine@monash.edu     
 */

#include <stdio.h>  
#include <unistd.h>   
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "yap.h"
#include "util.h"

long memallocd = 0;

double logadd(double V, double lp) {
  if ( lp>V ) {
    // swap so V is bigger
    double t = lp;
    lp = V;
    V = t;
  }
  return V + log(1.0+exp(lp-V));
}

unsigned samplet(float *fact, double facttot, unsigned dim, double U) {
  unsigned k;
  /*   got the relative proportions, so sample */
  if ( facttot<=0 ) {
    /*
     *    escape, just in case *everything* gets zeroed
     */
    k = U * dim;
  } else {
    facttot *= U;
    for (k=0; k<dim; k++) {
      facttot -= fact[k];
      if ( facttot< 0 )
        break;
    }
  }
  if ( k>=dim )
    k = dim - 1;
  return k;
}

uint16_t *u16vec(int n) 
{
  uint16_t *x = (uint16_t*)calloc(n,sizeof(uint16_t));
  memallocd += n*sizeof(uint16_t);
  if ( !x)
    yap_quit("cannot allocate uint16_t space of %ld\n", n*sizeof(uint16_t));
  return x;
}
uint32_t *u32vec(int n) 
{
  uint32_t *x = (uint32_t*)calloc(n,sizeof(uint32_t));
  memallocd += n*sizeof(uint32_t);
  if ( !x)
    yap_quit("cannot allocate uint32_t space of %ld\n", n*sizeof(uint32_t));
  return x;
}
  
double *dvec(int n) //
{
  double *x = (double*)calloc(n,sizeof(double));
  memallocd += n*sizeof(double);
  if ( !x)
    yap_quit("cannot allocate double space of %ld\n", n*sizeof(double));
  return x;
}
float *fvec(int n) //
{
  float *x = (float*)calloc(n,sizeof(float));
  memallocd += n*sizeof(float);
  if ( !x)
    yap_quit("cannot allocate float space of %ld\n", n*sizeof(float));
  return x;
}

#define BBLOCK 100000000
uint16_t **u16mat(int nr, int nc) //
{
  long int nrc = nr*(long)nc;
  uint16_t *tmp;
  uint16_t **x  = (uint16_t**)calloc(nr,sizeof(uint16_t*));
  memallocd += nrc*sizeof(uint16_t) + nr*sizeof(uint16_t*);
  unsigned r;
  /*
   *  nrc too big, can be 3Gb ina block,
   *  so we split into BBLOCK chunks
   */
  unsigned b, blockcnt, blocksize;
  if ( !x)
    yap_quit("cannot allocate uint16_t(mat) space of %ld\n", nr *sizeof(uint16_t*));
  blockcnt = (nrc-1)/BBLOCK + 1;
  if ( blockcnt>=nr ) 
    yap_quit("cannot allocate uint16_t(mat) space of %d*%d=%ld uint16s, too many blocks=%u\n", 
	     nr, nc, nrc, blockcnt);
  blocksize = (nr-1)/blockcnt + 1;
  // yap_message("u16mat(%d, %d): blocks=%d, size=%d\n", nr, nc, blockcnt, blocksize);
  for (b=0; b<blockcnt; b++) {
    /*
     *    allocate a chunk of tb columns
     */
    int tb = blocksize;
    if ( b==blockcnt-1 )
      tb = nr - b*blocksize;
    tmp = (uint16_t*) calloc(tb*(long)nc,sizeof(uint16_t));
    if ( !x)
      yap_quit("cannot allocate uint16_t(mat) space of %ld, total=%ld\n",  tb*nc*sizeof(uint16_t*),memallocd);
    for (r = 0; r < tb; r++) x[(b*blocksize)+r] = tmp + nc*r;
  }
  return x;
}

/*
 *  needs its own version of free
 */
void u16mat_free(uint16_t **mat, int nr, int nc) {
  unsigned b, blockcnt, blocksize;
  blockcnt = (nr*(long)nc-1)/BBLOCK + 1;
  blocksize = (nr-1)/blockcnt + 1;
  for (b=0; b<blockcnt; b++) 
    free(mat[b*blocksize]);
  free(mat);
}


uint32_t **u32mat(int nr, int nc) //
{
  long int nrc = nr*(long)nc;
  uint32_t *tmp = (uint32_t*) calloc(nrc,sizeof(uint32_t));
  uint32_t **x  = (uint32_t**)calloc(nr,sizeof(uint32_t*));
  unsigned r;
  memallocd += nrc*sizeof(uint32_t) + (nr*sizeof(uint32_t*));
  if ( !tmp || !x)
    yap_quit("cannot allocate uint32_t(mat) space of %ld\n", nrc *sizeof(uint32_t));
  for (r = 0; r < nr; r++) x[r] = tmp + nc*r;
  return x;
}

uint32_t ***u32tri(int ns, int nr, int nc) //
{
  long int nsrc = nr*nc*(long)ns;
  long int nsr = ns*(long)nr;
  uint32_t ***x  = (uint32_t***)calloc(ns,sizeof(uint32_t**));
  unsigned r, s;
  assert(nsrc>0);
  memallocd += nsrc*sizeof(uint32_t) + (nsr*sizeof(uint32_t*))
    + (ns*sizeof(uint32_t**));
  if ( !x)
    yap_quit("cannot allocate uint32_t-tri space of %ld\n", nsrc*sizeof(uint32_t));
  x[0] = (uint32_t**)calloc(nsr,sizeof(uint32_t*));
  if ( !x[0] )
    yap_quit("cannot allocate uint32_t-tri space of %ld\n", nsrc*sizeof(uint32_t));
  x[0][0] = (uint32_t*)calloc(nsrc,sizeof(uint32_t));
  if ( !x[0][0] )
    yap_quit("cannot allocate uint32_t-tri space of %ld\n", nsrc*sizeof(uint32_t));
  for (r = 1; r < nr; r++) 
    x[0][r] = x[0][r-1] + nc;
  for (s = 1; s < ns; s++) {
    x[s] = x[s-1] + nr;
    x[s][0] = x[s-1][0] + nr*nc;
    for (r = 1; r < nr; r++) 
      x[s][r] = x[s][r-1] + nc;
  }
  return x;
}

double **dmat(int nr, int nc) //
{
  long int nrc = nr*(long)nc;
  double *tmp = (double*) calloc(nrc,sizeof(double));
  double **x  = (double**)calloc(nr,sizeof(double*));
  int r;
  memallocd += nrc*sizeof(double) + (nr*sizeof(double*));
  if ( !tmp || !x)
    yap_quit("cannot allocate double(mat) space of %ld\n", nrc *sizeof(double));
  for (r = 0; r < nr; r++) x[r] = tmp + nc*r;
  return x;
}

float **fmat(int nr, int nc) //
{
  long int nrc = nr*(long)nc;
  float *tmp = (float*) calloc(nrc,sizeof(float));
  float **x  = (float**)calloc(nr,sizeof(float*));
  int r;
  memallocd += nrc*sizeof(float) + (nr*sizeof(float*));
  if ( !tmp || !x)
    yap_quit("cannot allocate float(mat) space of %ld\n", nrc *sizeof(float));
  for (r = 0; r < nr; r++) x[r] = tmp + nc*r;
  return x;
}

sparse_vec *smat(int nr) //
{
  sparse_vec *x = (sparse_vec *)calloc(nr,sizeof(sparse_vec));
  memallocd += nr*sizeof(sparse_vec);
  if ( !x)
    yap_quit("cannot allocate sparse_vec space of %ld\n", nr *sizeof(sparse_vec));
  return x;
}

unsigned countntot(char *fname) //
{
  unsigned ntot=0;
  char buf[BUFSIZ];
  FILE *fp = fopen(fname ,"r"); 
  if ( !fp ) 
    yap_sysquit( "Cannot open file '%s' for read\n", fname);
  while (fgets(buf, BUFSIZ, fp)) ntot++;
  fclose(fp);
  assert(ntot>0);
  return ntot;
}

static int countnnz16(int nr, int nc, uint16_t **x, uint16_t cutoff) //
{
  int i, j, nnz=0;
  for (i = 0; i < nr; i++) 
    for (j = 0; j < nc; j++) 
      if ( x[i][j]>cutoff ) nnz++;
  return nnz;
}
static int countnnz32tri(int ns, int nr, int nc, uint32_t ***x, uint32_t cutoff) //
{
  int i, j, k, nnz=0;
  for (k = 0; k < ns; k++) 
    for (i = 0; i < nr; i++) 
      for (j = 0; j < nc; j++) 
	if ( x[k][i][j]>cutoff ) nnz++;
  return nnz;
}
static int countnnz32(int nr, int nc, uint32_t **x, uint32_t cutoff) //
{
  int i, j, nnz=0;
  for (i = 0; i < nr; i++) 
    for (j = 0; j < nc; j++) 
      if ( x[i][j]>cutoff ) nnz++;
  return nnz;
}

void read_didwid(char *wfile, char *dfile, int N, uint32_t *d, uint32_t *w) //
{
  int i,wt,dt;
  FILE *fp1 = fopen(wfile ,"r"); 
  FILE *fp2 = fopen(dfile ,"r"); 
  if ( !fp1 ) 
    yap_sysquit( "Cannot open file '%s' for read\n", wfile);
  if ( !fp2 ) 
    yap_sysquit( "Cannot open file '%s' for read\n", dfile);
  for (i = 0; i < N; i++) {
    if (fscanf(fp1,"%d",&wt)) w[i]=wt-1;
    if (fscanf(fp2,"%d",&dt)) d[i]=dt-1;
  }
  if ( ferror(fp1) )
    yap_sysquit("Error on reading file '%s' ", wfile);
  if ( ferror(fp2) )
    yap_sysquit("Error on reading file '%s' ", dfile);
   fclose(fp1);
  fclose(fp2);
}

void read_u16vec(char *dfile, int N, uint16_t *d) //
{
  int i;
  FILE *fp = fopen(dfile ,"r"); 
  if ( !fp ) 
    yap_sysquit( "Cannot open file '%s' for read\n", dfile);
  for (i = 0; i < N; i++) {
    unsigned u;
    if ( !fscanf(fp," %u", &u) ) {
      yap_sysquit( "Cannot read from '%s' position %d\n", dfile, i);
    }
    d[i] = u;
  }
  if ( ferror(fp) )
    yap_sysquit("Error on reading file '%s' ", dfile);
  fclose(fp);
}

/*
 *   assumes is zero filled already
 */
void read_fmat(char *dfile, int N, int C, float **f) //
{
  float val;
  int n, c;
  FILE *fp = fopen(dfile ,"r"); 
  if ( !fp ) 
    yap_sysquit( "Cannot open file '%s' for read\n", dfile);
  while ( fscanf(fp,"%d %d %g",&n, &c, &val)==3 ) {
    if ( n>=0 && c>=0 && n<N && c<C )
      f[n][c] = val;
  }
  if ( ferror(fp) )
    yap_sysquit("Error on reading file '%s' ", dfile);
  fclose(fp);
}

void read_dvec(char *dfile, int N, double *d) //
{
  int i;
  FILE *fp = fopen(dfile ,"r"); 
  if ( !fp ) 
    yap_sysquit( "Cannot open file '%s' for read\n", dfile);
  for (i = 0; i < N; i++) {
    if ( !fscanf(fp," %lg",&d[i]) ) {
      yap_sysquit( "Cannot read from '%s' position %d\n", dfile, i);
    }
  }
  if ( ferror(fp) )
    yap_sysquit("Error on reading file '%s' ", dfile);
  fclose(fp);
}

void write_fvec(char *dfile, int N, float *d) //
{
  int i;
  FILE *fp = fopen(dfile ,"w");
  if ( !fp ) 
    yap_sysquit( "Cannot open file '%s' for write\n", dfile);
  for (i = 0; i < N; i++) {
    if ( !fprintf(fp,"%g\n",d[i]) ) {
      yap_sysquit( "Cannot write to '%s' position %d\n", dfile, i);
    }
  }
  if ( ferror(fp) )
    yap_sysquit("Error on writing file '%s' ", dfile);
  fclose(fp);
}
void write_dvec(char *dfile, int N, double *d) //
{
  int i;
  FILE *fp = fopen(dfile ,"w");
  if ( !fp ) 
    yap_sysquit( "Cannot open file '%s' for write\n", dfile);
  for (i = 0; i < N; i++) {
    if ( !fprintf(fp,"%lg\n",d[i]) ) {
      yap_sysquit( "Cannot write to '%s' position %d\n", dfile, i);
    }
  }
  if ( ferror(fp) )
    yap_sysquit("Error on writing file '%s' ", dfile);
  fclose(fp);
}

void write_dmat(char *dfile, int N, int C, double **d) //
{
  int i, j;
  FILE *fp = fopen(dfile ,"w");
  if ( !fp ) 
    yap_sysquit( "Cannot open file '%s' for write\n", dfile);
  for (i = 0; i < N; i++) {
    for (j = 0; j < C; j++) {
      if ( !fprintf(fp,"%d %d %lg\n",i, j, d[i][j]) ) {
	yap_sysquit( "Cannot write to '%s' position %d,\n", dfile, i, j);
      }
    }
  }
  if ( ferror(fp) )
    yap_sysquit("Error on writing file '%s' ", dfile);
  fclose(fp);
}

void write_fmat(char *dfile, int N, int C, float **f) //
{
  int i, j;
  FILE *fp = fopen(dfile ,"w");
  if ( !fp ) 
    yap_sysquit( "Cannot open file '%s' for write\n", dfile);
  for (i = 0; i < N; i++) {
    for (j = 0; j < C; j++) {
      if ( !fprintf(fp,"%d %d %g\n",i, j, f[i][j]) ) {
	yap_sysquit( "Cannot write to '%s' position %d,\n", dfile, i, j);
      }
    }
  }
  if ( ferror(fp) )
    yap_sysquit("Error on writing file '%s' ", dfile);
  fclose(fp);
}

void write_u16sparseco(int nr, int nc, uint16_t **x, char *fname,
		       uint16_t cutoff) //
{
  FILE *fp = fopen(fname,"w");
  int i, j;
  if ( !fp ) 
    yap_sysquit( "Cannot open file '%s' for write\n", fname);
  fprintf(fp, "%d\n", nr);
  fprintf(fp, "%d\n", nc);
  fprintf(fp, "%d\n", countnnz16(nr,nc,x,cutoff));
  for (i = 0; i < nr; i++) 
    for (j = 0; j < nc; j++) 
      if ( x[i][j]>cutoff ) 
	fprintf(fp, "%d %d %u\n", i, j, (unsigned)x[i][j]);
  if ( ferror(fp) )
    yap_sysquit("Error on writing file '%s' ", fname);
  fclose(fp);
}
void write_u16sparse(int nr, int nc, uint16_t **x, char *fname) {
  write_u16sparseco(nr, nc, x, fname, 0);
}

/*
 *  ignores entries not in sparse format,
 *  so make sure is zerod before use
 */
void read_u16sparse(int nr, int nc, uint16_t **x, char *fname) //
{
  FILE *fp = fopen(fname,"r");
  int i, u, v;
  int Nin;
  int sparse=0;
  if ( !fp )
    yap_sysquit("Cannot open file '%s' for read\n", fname);
  if ( fscanf(fp, "%d", &Nin)!=1 || Nin!=nr )
    yap_quit("Number rows wrong for sparse matrix in '%s', should be %d\n",
	     fname, nr);
  if ( fscanf(fp, "%d", &Nin)!=1 || Nin!=nc )
    yap_quit("Number columns wrong for sparse matrix in '%s', should be %d\n",
	     fname, nc);
  if ( fscanf(fp, "%d", &sparse)!=1 )
    yap_quit("Cannot read sparsity in matrix in '%s'\n", fname);
  for (i = 0; i < sparse; i++) {
    if ( fscanf(fp, "%d %d %u", &u, &v, &Nin) !=3 )
      yap_quit("Cannot read line %d in matrix in '%s'\n", i, fname);
    x[u][v] = Nin;
  }			       
  if ( ferror(fp) )
    yap_sysquit("Error on reading file '%s' ", fname);
  fclose(fp);
}
void read_u32sparse(int nr, int nc, uint32_t **x, char *fname) //
{
  FILE *fp = fopen(fname,"r");
  int i, u, v;
  int Nin;
  int sparse=0;
  if ( !fp )
    yap_sysquit("Cannot open file '%s' for read\n", fname);
  if ( fscanf(fp, "%d", &Nin)!=1 || Nin!=nr )
    yap_quit("Number rows wrong for sparse matrix in '%s', should be %d\n",
	     fname, nr);
  if ( fscanf(fp, "%d", &Nin)!=1 || Nin!=nc )
    yap_quit("Number columns wrong for sparse matrix in '%s', should be %d\n",
	     fname, nc);
  if ( fscanf(fp, "%d", &sparse)!=1 )
    yap_quit("Cannot read sparsity in matrix in '%s'\n", fname);
  for (i = 0; i < sparse; i++) {
    if ( fscanf(fp, "%d %d %u", &u, &v, &Nin) !=3 )
      yap_quit("Cannot read line %d in matrix in '%s'\n", i, fname);
    x[u][v] = Nin;
  }
  if ( ferror(fp) )
    yap_sysquit("Error on reading file '%s' ", fname);
  fclose(fp);
}

void read_u32sparsetri(int ns, int nr, int nc, uint32_t ***x, char *fname) 
{
  FILE *fp = fopen(fname,"r");
  int i, u, v, w;
  int Nin;
  int sparse=0;
  if ( !fp )
    yap_sysquit("Cannot open file '%s' for read\n", fname);
  if ( fscanf(fp, "%d", &Nin)!=1 || Nin!=ns )
    yap_quit("Number blocks wrong for sparse matrix in '%s', should be %d\n",
	     fname, ns);
  if ( fscanf(fp, "%d", &Nin)!=1 || Nin!=nr )
    yap_quit("Number rows wrong for sparse matrix in '%s', should be %d\n",
	     fname, nr);
  if ( fscanf(fp, "%d", &Nin)!=1 || Nin!=nc )
    yap_quit("Number columns wrong for sparse matrix in '%s', should be %d\n",
	     fname, nc);
  if ( fscanf(fp, "%d", &sparse)!=1 )
    yap_quit("Cannot read sparsity in matrix in '%s'\n", fname);
  for (i = 0; i < sparse; i++) {
    if ( fscanf(fp, "%d %d %d %u", &u, &v, &w, &Nin) !=4 )
      yap_quit("Cannot read line %d in matrix in '%s'\n", i, fname);
    x[u][v][w] = Nin;
  }			       
  if ( ferror(fp) )
    yap_sysquit("Error on reading file '%s' ", fname);
  fclose(fp);
}

void write_u32sparseco(int nr, int nc, uint32_t **x, char *fname,
		       uint32_t cutoff ) //
{
  FILE *fp = fopen(fname,"w");
  int i, j;
  if ( !fp )
    yap_sysquit("Cannot open file '%s' for write\n", fname);
  fprintf(fp, "%d\n", nr);
  fprintf(fp, "%d\n", nc);
  fprintf(fp, "%d\n", countnnz32(nr,nc,x,cutoff));
  for (i = 0; i < nr; i++) 
    for (j = 0; j < nc; j++) 
      if ( x[i][j]>cutoff ) 
	fprintf(fp, "%d %d %u\n", i, j, (unsigned)x[i][j]);
  if ( ferror(fp) )
    yap_sysquit("Error on writing file '%s' ", fname);
  fclose(fp);
}
void write_u32sparse(int nr, int nc, uint32_t **x, char *fname) {
  write_u32sparseco(nr, nc, x, fname, 0);
}

void write_u32sparsetri(int ns, int nr, int nc, uint32_t ***x, char *fname,
			uint32_t cutoff) //
{
  FILE *fp = fopen(fname,"w");
  int i, j, k;
  if ( !fp )
    yap_sysquit("Cannot open file '%s' for write\n", fname);
  fprintf(fp, "%d\n", ns);
  fprintf(fp, "%d\n", nr);
  fprintf(fp, "%d\n", nc);
  fprintf(fp, "%d\n", countnnz32tri(ns,nr,nc,x,cutoff));
  for (k = 0; k < ns; k++) 
    for (i = 0; i < nr; i++) 
      for (j = 0; j < nc; j++) 
	if ( x[k][i][j]>cutoff ) 
	  fprintf(fp, "%d %d %d %u\n", k, i, j, (unsigned)x[k][i][j]);
  if ( ferror(fp) )
    yap_sysquit("Error on writing file '%s' ", fname);
  fclose(fp);
}

void write_u16vec (int n, uint16_t *x, char *fname) //
{
  FILE *fp = fopen(fname,"w");
  int i;
  if ( !fp )
    yap_sysquit("Cannot open file '%s' for write\n", fname);
  for (i = 0; i < n; i++)  
    fprintf(fp, "%d\n", x[i] );
  if ( ferror(fp) )
    yap_sysquit("Error on writing file '%s' ", fname);
  fclose(fp);
}

void chksum (int n, int T, double **x, double *sumx) //
{
  int i, t;
  double sum;
  for (t = 0; t < T; t++) {
    sum = 0.0;
    for (i = 0; i < n; i++) sum += x[i][t];
    assert(fabs(sum-sumx[t])<1e-6);
  }
}

double  dmax(int n, double *x) //
{
  int i;
  double xmax=x[0];
  for (i = 0; i < n; i++) xmax = MAX(xmax,x[i]);
  return xmax;
}
double  dmin(int n, double *x) //
{
  int i;
  double xmax=x[0];
  for (i = 0; i < n; i++) xmax = MIN(xmax,x[i]);
  return xmax;
}

uint32_t u32max(int n, uint32_t *x) //
{
  int i;
  uint32_t xmax=x[0];
  for (i = 0; i < n; i++) xmax = MAX(xmax,x[i]);
  return xmax;
}

void vec_atod(double *vec, int N, char *vals) {
  char *val;
  int n = 0;
  while ( n<N && (val=strsep(&vals," ")) ) {
    sscanf(val, " %lg", &vec[n]);
    n++;
  }
  for ( ; n<N; n++) vec[n] = 0;
}

/*
 *  pick up line from file "stem.ext" starting with "par" and an "=",
 *  and return stuff after "=" but only first len chars
 */
static char *readextpar(char *stem, char *ext, char *par, char *buf, int len) {
  char *file = yap_makename(stem,ext);
  FILE *fp = fopen(file,"r");
  int parlen = strlen(par);
  int buflen = len+parlen+50;
  char *linebuf = malloc(buflen+2);
  if ( !fp ) 
    yap_sysquit("Cannot open parameter file '%s' ", file);
  if ( !linebuf )
    yap_quit("Out of memory in readpar(%s)\n", par);
  buf[0] = 0;
  linebuf[0] = 0;
  while ( fgets(&linebuf[0],buflen,fp) ) {
    if ( strstr(&linebuf[0],par)==&linebuf[0] && 
	 (linebuf[parlen]==' ' || linebuf[parlen]=='=') ) {
      char *ret = strstr(&linebuf[0],"=");
       if ( !ret ) {
	yap_quit("Badly formatted parameter file '%s': %s\n",
		 file, &linebuf[0]);
      }
      ret++;
      // yap_message("ret = '%s', length %d\n", ret, strlen(ret));
      {
        int i;
        i = strlen(ret);
        if ( i>=len ) {
          i = len-1;
          ret[i] = 0;
          i--;
        }
        for ( ; i>=0; i--)
          buf[i] = ret[i];
        /*  for some weird reason, this fails occasionally!! */
        // strncpy(buf, ret, len);
      }
      strncpy(buf, ret, len);
      break;
    }
  }
  if ( ferror(fp) )
    yap_sysquit("Error on parameter file '%s' ", file);
  fclose(fp);
  free(file);
  free(linebuf);
  if ( buf[0] ) 
    return &buf[0];
  else
    return NULL;
}

char *readpar(char *stem, char *par, char *buf, int len) {
  return readextpar(stem,".par",par,buf,len);
}
char *readsrcpar(char *stem, char *par, char *buf, int len) {
  return readextpar(stem,".srcpar",par,buf,len);
}
