/*
 * Lower level routines on statistics
 * Copyright (C) 2010-2014 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@monash.edu)
 *     
 */

#ifndef __CHANGE_H
#define __CHANGE_H

#include "hca.h"

/*
 *  modifying indicators for topic or word during Gibbs sampling
 */
int tableidtopic(int d, int t, int tot);
void fix_tableidtopic(int d, int t);
void fix_tableidword(int w, int t);
void unfix_tableidtopic(int d, int t);
void unfix_tableidword(int w, int t);

/*
 *  adjusting stats for entire doc during testing
 */
int remove_doc(int did, enum GibbsType fix);
int add_doc(int did, enum GibbsType fix);
void zero_doc(int did);

/*
 *  computing stats
 */
uint16_t comp_Td(int did);
int nonzero_Ndt(int t);
int nonzero_Nwt(int t);

#endif
