/*
 * Atomic operations wrapper
 * Copyright (C) 2014 Wray Buntine and Swapnil Mishra
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@monash.edu)
 *         Swapnil Mishra
 *
 *
 *  Note return values for the _fetch() routines currently *NOT* used
 *  so can do  _atomic_add_fetch() or _atomic_fetch_add() 
 */
#ifndef __ATOMIC_H
#define __ATOMIC_H

extern long atomic_zero;

/*
 *   define to make things non-atomic
 */
// #define NONATOMIC

#ifdef NONATOMIC
#define atomic_incr(inttype) ++(inttype)
/*
 *   if its zero, incr. and return true, else do nothing and return false
 */
#define atomic_incr_zero(inttype) ((inttype==0)?(inttype)++:0)
#define atomic_decr(inttype) --(inttype)
#define atomic_add(inttype,val) (inttype += val)
#define atomic_sub(inttype,val) (inttype -= val)
#else
/* Test for GCC == 4.8.2 */
#if (__GNUC__==4 && __GNUC_MINOR__==8 && \
     ( __GNUC_PATCHLEVEL__<=2 && __GNUC_PATCHLEVEL__>=1) )
#define atomic_incr_zero(inttype)  ((inttype==0)?atomic_incr(inttype):0)
#define atomic_incr(inttype) __atomic_add_fetch(&(inttype),1, __ATOMIC_RELAXED)
#define atomic_decr(inttype) __atomic_sub_fetch(&(inttype),1, __ATOMIC_RELAXED)
#define atomic_add(inttype,val) __atomic_add_fetch(&(inttype),val, __ATOMIC_RELAXED)
#define atomic_sub(inttype,val) __atomic_sub_fetch(&(inttype),val, __ATOMIC_RELAXED)
#else
#if (__GNUC__==4 && __GNUC_MINOR__==1 &&__GNUC_PATCHLEVEL__==2  )
#define atomic_incr_zero(inttype) ???
#define atomic_incr(inttype) __sync_add_and_fetch(&(inttype),1)
#define atomic_decr(inttype) __sync_sub_and_fetch(&(inttype),1)
#define atomic_add(inttype,val) __sync_add_and_fetch(&(inttype),val)
#define atomic_sub(inttype,val) __sync_sub_and_fetch(&(inttype),val)
#else
/*
 *  leave undefined to force non compile
 */
#define atomic_incr_zero(inttype) ???
#define atomic_incr(inttype) ???
#define atomic_decr(inttype) ???
#define atomic_add(inttype,val) ???
#define atomic_sub(inttype,val) ???
#endif
#endif
#endif

#endif
