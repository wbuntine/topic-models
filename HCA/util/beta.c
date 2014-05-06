/*
 * 
 * Copyright (C) 2011 Wray Buntine
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (Wray.Buntine@nicta.com.au)
 *
 */

#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "srng.h"
#include "yap.h"
extern rngp_t rngp;

double ran_beta(double b, int n) {
  int try = 50;
  double q = 0;
  while ( --try>=0 && q<=0 ) {
    /*
     *   occasionally when b is small it gets underflow
     *   so we have to redo
     */
    q = rng_beta(rngp, b, n);
  }
  if ( q<=0 ) {
    yap_message("q in ran_beta(b=%f,n=%d) went zero\n", b, n);
    return 0;
  }
  return q;
}

