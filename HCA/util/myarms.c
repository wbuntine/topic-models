/*
 * simpler interface to arms.h
 * with checking of bounds
 */


#include "stdlib.h"
#include "math.h"
#include "arms.h"
#include "yap.h"

#ifdef isfinite
#define ISFINITE(x) isfinite(x)
#else
#define ISFINITE(x) finite(x)
#endif

int myarms_evals;
double myarms_last;
extern int verbose;

static int myarms_simple(int ninit, double *xl, double *xr,
			 double (*myfunc)(double x, void *mydata), void *mydata,
			 int dometrop, double *xprev, double *xsamp, int nsamp)

/* adaptive rejection metropolis sampling - simplified argument list */
/* ninit        : number of starting values to be used */
/* *xl          : left bound */
/* *xr          : right bound */
/* *myfunc      : function to evaluate log-density */
/* *mydata      : data required by *myfunc */
/* dometrop     : whether metropolis step is required */
/* *xprev       : current value from markov chain */
/* *xsamp       : to store sampled value */
/* nsamp        : number of samples */

{
  double xinit[ninit+1], convex=1.0, qcent, xcent;
  int err, i, npoint=100, ncent=0, neval;

  /* set up starting values */
  for(i=0; i<ninit; i++){
    xinit[i] = *xl + (i + 0.1) * (*xr - *xl)/(ninit - 1.0 + 0.2);
  }
  /*  now insert *xprev */
  if ( *xprev>*xl && *xprev<*xr ) {
    /*  now insert *xprev */
    for(i=0; i<ninit && xinit[i]<*xprev; i++) ;
    if ( i>=ninit )
      xinit[ninit] = *xprev;
    else {
      int savei = i;
      for(i=ninit; i>savei; i--) {
	xinit[i] = xinit[i-1];
      }
      xinit[savei] = *xprev;
    }
  }
  err = arms(xinit,ninit,xl,xr,myfunc,mydata,&convex,npoint,
	     dometrop,xprev,xsamp,
             nsamp,&qcent,&xcent,ncent,&neval);
  return err;
}

#define NSAMP 10
int myarmsMH(double xl, double xr,
	     double (*myfunc)(double x, void *mydata), 
	     void *mydata, double *xval, char *label,
	     int doMH) {
  double result = *xval;
  double *resvec = NULL;
  double startval = *xval;
  int errcode;
  if ( fabs(*xval-xr)/(xr)<0.00001 ) {
    *xval = xr * 0.999 + xl * 0.001;
  }
  if ( fabs(*xval-xl)/xl<0.00001 ) {
    *xval = xl * 0.999 + xr * 0.001;
  }
  myarms_evals = 0;
  if ( doMH ) {
    resvec = malloc(sizeof(resvec[0])*NSAMP);
    errcode = myarms_simple(6, &xl, &xr, myfunc, mydata, doMH, xval, 
			  resvec, NSAMP);
    result = resvec[NSAMP-1];
    free(resvec);
  } else 
    errcode = myarms_simple(6, &xl, &xr, myfunc, mydata, 0, xval, &result, 1);
  /*
   *  1007, 1003 is out of bounds
   */
  if ( errcode && errcode!=1007 && errcode!=1003
       && (errcode!=2000 || startval!=result  ) ) {
    yap_quit("   myarmsMH(%lf,%s)->%d = %lf,%lf%s->%lf, w %d calls, quitting\n", 
	     startval, label, errcode,
	     myarms_last, result, (!ISFINITE(result))?"(inf)":"",
	     *xval, myarms_evals);
  }
  if ( errcode==1007 || errcode==1003 ) {
	  yap_message("   error myarmsMH(%lf,%s)->%d = %lf,%lf%s->%lf, w %d calls, quitting\n",
             startval, label, errcode,
             myarms_last, result, (!ISFINITE(result))?"(inf)":"",
             *xval, myarms_evals);
          /*  hit bounds, for safety use start val */
	  result = startval;
  } 
  /*
   *    note, sometimes the value is returned
   *    unchanged .... seems to be when the 
   *    posterior is really focussed
   */
  if ( verbose>1 || !ISFINITE(result) || startval==result )
    yap_message("   myarmsMH(%s) = %lf%s<-%lf, w %d calls %s\n", 
		label,
		result, (!ISFINITE(result))?"(inf)":"",
		*xval, myarms_evals,
		(startval==result)?"UNCHANGED":"");
  if ( ISFINITE(result) ) {
    if ( result<xl )
      result = xl;
    else if ( result>xr )
      result = xr;
    *xval = result;
  }
  if ( !ISFINITE(result) ) {
    return 1;
  }
  return 0;
}

int myarms(double xl, double xr,
	   double (*myfunc)(double x, void *mydata), 
	   void *mydata, double *xval, char *label) {
  return myarmsMH(xl, xr, myfunc, mydata, xval, label, 0);
}

