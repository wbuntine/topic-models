/*
 * PMI computation from top topics previously computed
 * Copyright (C) 2013 Wray Buntine
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine
 */
#ifndef __PMI_H
#define __PMI_H

/****************************************
 *  reading of topics file
 */
void ttop_open(char *topfile);
void ttop_close();
/*
 *    get next (topic,epoch) pair
 *    (T,E) are input maximum topics and epochs
 *    return zero if none
 *    otherwise return topic and epoch in (*k,*e)
 */
int ttop_next(int T, int E, int *k, int *e);
/*
 *    get next word
 *    return zero if none
 *    returned in *j
 *    W is bound for word index
 *    n_word is the position for this list ... used for error reporting
 */
int ttop_word(int W, int n_word, unsigned int *j);
/*
 *  free buffer
 */
void ttop_eol();


/*
 *    if tpmi==NULL
 *         print out PMI for topics computed on topk words
 *    else
 *         store in tpmi
 *    tpmi[e*(ddN.T+1)+k] = PMI for topic k in epoch e
 *    tpmi[e*(ddN.T+1)+ddN.T] = mean (by topic propportion) PMI for epoch 
 *
 *    return average over epochs
 */
double report_pmi(char *topfile,   /* name of topics file */
		  char *pmifile,  /* name of PMI file */
                  char *toppmifile, /* name of topics+pmi file to write */
		  int T,          /* total topics */
		  int W,          /* total words */
		  int E,          /*  number of epochs */
		  int topk,
		  double *tp,
                  float *tpmi);

#endif
