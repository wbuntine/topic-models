/*
 *   YAP is for reporting errors.
 *   By default, reports will be time-stamped and sent to stderr.
 *   Copyright (C) 2002-2006 Wray Buntine
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (buntine@hiit.fi)
 *
 */

#ifndef _YAP_H_
#define _YAP_H_
#include <stdarg.h>

/*
 *    few general things here since yap.h gets included a lot
 */
#ifdef __cplusplus
#  define EXTERN        extern "C"
#else
#  define EXTERN        extern
#endif

#if (__GNUC__ >2 || __GNUC_MINOR__ >=7) && !defined(UNUSED)
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#if defined(__STDC__) || defined(__cplusplus) || defined(_MSC_VER)
#define __STRING(x)     #x
#else   /* !(__STDC__ || __cplusplus) */
#define __STRING(x)     "x"
#error __STRING is not supported!
#endif

#define yap_free(ptr)  if (!ptr) yap_message("Freeing NULL at %s:%d\n",__FILE__, __LINE__); else free(ptr);

#ifdef NDEBUG
#define yap_infinite(dbl)  
#else
#define yap_infinite(dbl)  if (!finite(dbl)) yap_message("Var '%s' got infinite at %s:%d\n",__STRING(dbl),__FILE__, __LINE__);
#endif


#define yap_free2(ptr)  

#define yap_realloc(ptr,oldsize,size) realloc(ptr,size) 

/*
 *   A realloc that uses malloc, but ONLY works if increasing in size
 */
extern void *yap_fullrealloc(void *ptr, size_t oldsize, size_t size);

/*
 *   Use syslog(,LOG_LOCAL2) to report errors. 
 */
extern void yap_log(char *name);
/*
 *  some details to STDERR
 */
extern void yap_about();
/* 
 *  Use write()s to this fd to report errors. 
 *  Note writes() means no buffering.
 *  The file will be opened in append mode.
 */
extern void yap_file(char *name);
/* 
 *  Use this to stamp messages with time/pid etc. 
 *  yap_stamp_in() is your routine to create the string
 *  and yap assumes nothing about memory for it, just
 *  calls and copies.
 */
extern void yap_setstamp(char *(*yap_stamp_in)(void));
/*
 *  the default for the above
 */
extern char *yap_stamp(void);
/*
 *  assumes error has been redirected somehow;
 *  writes out the PID and goes in background
 *  though doesn't change directory;
 *  set watchdog to put a watchdog restarter in place
 */
void yap_daemon(char *pidfilename, int watchdog);

/* Nonfatal error related to a system call. */
extern void yap_sysreport(const char *fmt, ...);
/* Nonfatal error unrelated to a system call. */
extern void yap_report(const char *fmt, ...);
/* Nonfatal message with no additional error/stamp. */
extern void yap_message(const char *fmt, ...);
extern void yap_va(const char *fmt, va_list ap);
/* Fatal error related to a system call. */
extern void yap_sysquit(const char *fmt, ...);
/* Fatal error unrelated to a system call. */
extern void yap_quit(const char *fmt, ...);
/* print UNIX calling command line using yap_report() */
extern void yap_commandline(int argc, char**argv);

extern char *yap_makename(char *stem, char *ext);
extern char *yap_makenamedot(char *stem, char *ext);
extern int  yap_fileexists(char *fname);
#endif
