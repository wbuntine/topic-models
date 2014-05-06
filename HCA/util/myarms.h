/*
 * simpler interface to arms.h
 * with checking of bounds
 */

extern int myarms_evals;
extern double myarms_last;

int myarms(double xl, double xr,
	   double (*myfunc)(double x, void *mydata), 
	   void *mydata, double *xval, char *label);
int myarmsMH(double xl, double xr,
	     double (*myfunc)(double x, void *mydata), 
	     void *mydata, double *xval, char *label,
	     int doMH);
