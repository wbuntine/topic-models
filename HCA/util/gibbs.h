/*
 * Gibbs help for LDA style latent variable modelling
 * Copyright (C) 2014 Wray Buntine
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@monash.edu)
 *
 *
 */
#ifndef __GIBBS_H
#define __GIBBS_H

/*
 *  fix==GibbsHold for document completion (or word hold-out) testing
 *       as it stops some words affecting topic
 *  fix==GibbsNone for regular training
 */
enum GibbsType { GibbsNone, GibbsHold };

/*
 *  topics are stored in uint16_t's
 *  and keep a Boolean called 'r' in high bit to indicate
 *  if the word is bursty
 */
#define Z_r(z) ((z&32768U)?1:0)
#define Z_t(z) (z&32767U)
#define Z_issetr(z) (z&32768U)
#define Z_setr(z) z |= 32768U
#define Z_unsetr(z) z &= 32767U
#define Z_sett(z,val) {z &= 32768U; z |= (uint16_t)(val);}

#endif
