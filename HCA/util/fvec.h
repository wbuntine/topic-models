/*
 * Utilities
 * Copyright (C) 2014 Wray Buntine
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@monash.edu     
 */

void fv_copy(float *v1, float *v2, int N) ;
/*
 *  Hellinger distance fr prob vecs
 *    \sum_i (sqrt(p[i])-sqrt(q[i]))^2 
 *     = 2 - 2\sum_i sqrt(p[i]q[i])
 */
double fv_helldistunif(float *vp, int N) ;
double fv_helldist(float *vp, float *vp2, int N) ;
double fv_entropy(float *vp, int N) ;
double fv_kl(float *vp1, float *vp2, int N);
double fv_expprob(float *vp, int N);
double fv_avestrlen(float *vp, char **str, int N) ;
double fv_bound(float *vp, int N, float alpha);
