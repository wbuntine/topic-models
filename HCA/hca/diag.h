/*
 * Various data for diagnostics
 * Copyright (C) 2013 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *     
 */

#ifndef __DIAG_H
#define __DIAG_H

/*
 *  diagnostic and parameter estimates kept by the training algorithm
 */
typedef struct D_diag_s {
  /*
   *     used to compute docXtopic probs .... is a memory hog
   *     separate for training and test (tprob is test version)
   */
  uint16_t    didtprob;
  uint16_t    didprob;
  uint16_t    doprob;
  float **prob;
  float **tprob;
  /*
   *    for estimation of topicXword \phi matrix
   */
  uint16_t    dophi;
  uint16_t    doalpha;
  int     phi_cnt;
  int     alpha_cnt;
  /*
   *   used to compute sparsity over time
   */
  uint32_t *iscodeword;
  uint32_t *words;
  int     n_words;
  uint16_t    docode;
  uint16_t    didcode;
  float ***code;
} D_diag_t;


#define G_isword(w) (ddG.iscodeword[(w)/32U] & (1U << (((unsigned)w)%32U)))

extern D_diag_t ddG;

void diag_alloc();

/*
 *  optionally save/update \phi matrix at end of cycles
 */
void phi_init(char *resstem);
void phi_update();
void phi_save();
void phi_free();
void phi_load(char *resstem);
double phi_entropy(int k);

/*
 *  optionally save/update \alpha matrix at end of cycles
 */
void alpha_init(char *resstem);
void alpha_update();
void alpha_save();
void alpha_free();
double alpha_entropy();

/*
 *    during Gibbs estimate proportion
 *    of topic in this select words
 */
void sparsemap_init(FILE *fp, int proc);
void sparsemap_null();
void sparsemap_free();
void sparsemap_report(char *resstem, double epsilon, int proc);
int sparsemap_word(uint32_t w);

/*
 *    during Gibbs estimate proportion
 *    of topics for all documents
 */
void tprob_init();
void tprob_null();
void tprob_free();
void tprob_report(char *resstem, double epsilon);
void prob_report(char *resstem, double epsilon);
void prob_load(char *resstem, char *suff, float **mat);

#endif
