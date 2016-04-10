/*
 * Various checking routines for debugging
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

#ifndef __CHECK_H
#define __CHECK_H

/*
 *  checks
 */
void check_Ndt(int d);
void check_Tw();
void check_Nwt(int w, int val);
void check_TWT();
void check_TDT();
void check_sparse();

#endif
