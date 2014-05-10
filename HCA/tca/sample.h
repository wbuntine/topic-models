#include <unistd.h>

#ifndef __SAMPLE_H
#define __SAMPLE_H

void sample_bm0(double *b);
void sample_bm1(double *b);
void sample_bp0(double *b);
void sample_bp1(double *b);
void sample_bt(double *b);
void sample_bb(double *b);
void sample_at(double *a);
void sample_am(double *a);
void sample_ap0(double *a);
void sample_ap1(double *a);
void sample_ab(double *a);


extern int verbose;
double likelihood();
void cache_update(char *par);
double ran_beta(double b, int n);

#endif
