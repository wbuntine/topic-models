/*
 * Various bounds and initialisations for DPs. PYPs and Dirs
 * Copyright (C) 2014 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@monash.au)
 *     
 */

#ifndef __DIRDM_H
#define __DIRDM_H

/*
 *    bounds for the parameters to various distributions
 */
#define DIR_MIN 0.0001
#define DIR_MAX 1.0
#define DIR_TOTAL_MAX 10000.0
#define PYP_DISC_MIN 0.01
#define PYP_DISC_MAX 0.98
#define PYP_CONC_MIN 0.001
#define PYP_CONC_MAX 10000.0

/*
 *    priors for PYP_CONC
 *    NB.  scale is multiplied by the dimension
 */
#define PYP_CONC_PSHAPE 1.1
#define PYP_CONC_PSCALE 1

#endif
