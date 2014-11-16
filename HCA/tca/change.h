/*
 * Changing data structures for statistics
 * Copyright (C) 2014 Jinjing Li and Wray Buntine
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Authors: Jinjing Li <jinjingli@gmail.com>
 *          Wray Buntine (wray.buntine@monash.edu)
 *          
 */

#ifndef __CHANGE_H
#define __CHANGE_H

#include "gibbs.h"

/*
 *  changing latent stats
 */
void fix_tableidtopic(int d, int t, int ind);
void unfix_tableidtopic(int d, int t, int ind);
/*
 *   these also return the smallest epoch e for which
 *   the phi_sum_cache[*][t][e] is invalidated
 */
int fix_tableidword(int e, int w, int t, int ind);
int unfix_tableidword(int e, int w, int t, int ind);

/*
 *  computing stats
 */
uint16_t comp_Td(int did);
int nonzero_n_dt(int t);
int nonzero_m_vte(int e, int t);

/*
 *  changing all stats for doc
 */
int remove_doc(int did, enum GibbsType fix);
int add_doc(int did, enum GibbsType fix);

#endif
