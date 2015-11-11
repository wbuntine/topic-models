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

/*
 *   define to make things non-atomic
 */
// #define NONATOMIC
#ifndef H_THREADS
#define NONATOMIC
#endif


#ifdef NONATOMIC
/* 
 *    *no* atomic ops
 */
#define atomic_incr(inttype) (++(inttype))
#define atomic_decr(inttype) (--(inttype))
/*
 *   if its val, incr./decr. and return true, else do nothing and return false
 */
#define atomic_incr_val(inttype,val) ((inttype==val)?((inttype)++,1):0)
#define atomic_decr_val(inttype,val) ((inttype==val)?((inttype)--,1):0)
#define atomic_add(inttype,val) (inttype += val)
#define atomic_sub(inttype,val) (inttype -= val)
#else
#if ( (__GNUC__==4 && __GNUC_MINOR__>=7) || __GNUC__>4 )
/* 
 *    Test for GCC == 4.[789].? 
 */
#define atomic_decr_val(inttype,val)  __atomic_compare_exchange_n(&(inttype),&val,(val-1),0,__ATOMIC_RELAXED,__ATOMIC_RELAXED)
#define atomic_incr_val(inttype,val)  __atomic_compare_exchange_n(&(inttype),&val,(val+1),0,__ATOMIC_RELAXED,__ATOMIC_RELAXED)
#define atomic_incr(inttype) __atomic_add_fetch(&(inttype),1, __ATOMIC_RELAXED)
#define atomic_decr(inttype) __atomic_sub_fetch(&(inttype),1, __ATOMIC_RELAXED)
#define atomic_add(inttype,val) __atomic_add_fetch(&(inttype),val, __ATOMIC_RELAXED)
#define atomic_sub(inttype,val) __atomic_sub_fetch(&(inttype),val, __ATOMIC_RELAXED)
#else
#if (__GNUC__==4 && (( __GNUC_MINOR__==1 &&__GNUC_PATCHLEVEL__==2) || ( __GNUC_MINOR__>=4 && __GNUC_MINOR__<=6 ) )  )
/* 
 *    Test for GCC == 4.1.2 or  GCC == 4.4.?
 */
#define atomic_decr_val(inttype,val) __sync_bool_compare_and_swap(&(inttype),val,(val-1))
#define atomic_incr_val(inttype,val) __sync_bool_compare_and_swap(&(inttype),val,(val+1))
#define atomic_incr(inttype) __sync_add_and_fetch(&(inttype),1)
#define atomic_decr(inttype) __sync_sub_and_fetch(&(inttype),1)
#define atomic_add(inttype,val) __sync_add_and_fetch(&(inttype),val)
#define atomic_sub(inttype,val) __sync_sub_and_fetch(&(inttype),val)
#else
/*
 *  leave undefined to force non compile
 */
#define atomic_incr_val(inttype,val) ???
#define atomic_decr_val(inttype,val) ???
#define atomic_incr(inttype) ???
#define atomic_decr(inttype) ???
#define atomic_add(inttype,val) ???
#define atomic_sub(inttype,val) ???
#endif
#endif
#endif

#endif
