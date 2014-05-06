#ifndef __SAMPLE_H
#define __SAMPLE_H

#include "cache.h"


void sample_b(double *bw);
void sample_a(double *aw);
void sample_a0(double *aw0);
void sample_b0(double *bw);
void sample_bw(double *bw);
void sample_aw(double *aw);
void sample_aw0(double *aw0);
void sample_bw0(double *bw);
void sample_alpha(double *alpha);
void sample_beta(double *mytbeta);
void sample_ad(double *ad);
void sample_adk(double *ad);
void sample_bd(double *bd);
void sample_bdk(double *bdk);

extern int verbose;
double likelihood();
double ran_beta(double b, int n);

#endif
