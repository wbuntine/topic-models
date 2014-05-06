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
#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <signal.h>

// *nix only headers for windows use yapwin.c instead
#ifndef HAVE_CYGWIN     /* Cygwin does not need this */
#include <libgen.h>
#endif
#include <unistd.h>
#include <syslog.h>

#include "yap.h"

/*  flags for yap_print() */
#define YAP_ERRNO 2
#define YAP_MSG 4

/*  forward decl. */
extern void yap_sysquit(const char *fmt, ...);

/*
 *  Error logging by default goes to standard error.
 */
static char progname[20] = "";
static int  use_syslog = 0;
static int  yap_fd = STDERR_FILENO;

/*
 *  user defined or default error stamp
 */
static char *(*use_yap_stamp)(void) = yap_stamp;


/*
 * Use syslog(,LOG_LOCAL2) to report errors.
 * name = program name, only first 19 chars used
 */
void yap_log(char *name)
{
  strncpy(progname,name,19);
  progname[19] = 0;
  openlog(progname, LOG_PID|LOG_NDELAY, LOG_LOCAL2);
  use_syslog = 1;
  if ( yap_fd > 0 ) {
    close(yap_fd);
    yap_fd = -1;
  }
}

void yap_about() {
  fprintf(stderr,"YAP:\nprogname=%s\n", progname);
  fprintf(stderr,"yap_fd=%d vs. %d\n", yap_fd, STDERR_FILENO);
  fprintf(stderr,"use_syslog=%d\n", use_syslog);
}

/*
 * Use write()s to this fd to report errors.
 * Safe to call if simply reopening.
 * name = name of error file
 */
void yap_file(char *name)
{
  int fd;
  /*  put error log into safe mode beforehand */
  if ( yap_fd>0 && yap_fd!=STDERR_FILENO) {
    close(yap_fd);
    yap_fd = STDERR_FILENO;
  }
  if ((fd = open(name, O_CREAT | O_WRONLY | O_APPEND, 0644)) < 0) {
    yap_sysquit("ERROR:  opening %s", name);
  }
  yap_fd = fd;
  yap_report("Error logging opened on '%s' file number %d\n",
	     name, fd);
  use_syslog = 0;
}

static pid_t yapchildpid = -1;
static char *yappiffile = NULL;

static void doonexit(int sig) {
  if ( yapchildpid>1 ) {
    if ( kill(yapchildpid,SIGTERM)<0 )
      yap_sysreport("Cannot kill child %d\n", yapchildpid);
    else
      yap_report("Killed child %d\n", yapchildpid);
  } 
  if ( yappiffile ) {
    unlink(yappiffile);
  }
  exit(0);
}

void yap_daemon(char *pidfilename, int watchdog) {
  FILE *dp;
  if ( yap_fd == STDERR_FILENO && !use_syslog )
    yap_message("No error logging set up\n");
  errno = 0;
  if ( !daemon(1,0) && errno )
    yap_sysquit("Cannot go in background as daemon\n");
  yappiffile = strdup(pidfilename);
  dp = fopen(pidfilename,"w");
  if ( !dp )
    yap_sysquit("Cannot write to '%s'\n",pidfilename);
  fprintf(dp, "%d\n", getpid());
  fclose(dp);
  if ( watchdog ) {
    int stat_loc;
    while ( (yapchildpid = fork()) ) {
      if ( signal(SIGHUP,doonexit)==SIG_ERR )
	yap_sysquit("Cannot set signal\n");
      if ( signal(SIGTERM,doonexit)==SIG_ERR )
	yap_sysquit("Cannot set signal\n");
      /*
       *  parent waits around to restart when needed
       *  NB.  child seems to inherit sig actions too ??
       */
      if ( !waitpid(yapchildpid,&stat_loc,0) )
	yap_sysquit("yap_daemon failed waitpid\n");
      yapchildpid = -1;
      /*
       *   if an OK termination, then quit
       */
      if ( !stat_loc ) {
	exit(0);
      }
      yap_report("Watchdog restarting\n");
      signal(SIGHUP,SIG_DFL);
      signal(SIGTERM,SIG_DFL);  
    }
  }
}

void yap_setstamp(char *(*yap_stamp_in)(void))
{
  if ( yap_stamp_in )
    use_yap_stamp = yap_stamp_in;
  else 
    use_yap_stamp = yap_stamp;
}


/*
 * Return static string containing current time.
 */
char *yap_stamp(void)
{
  static char   timestring[28];
  struct tm     *tmp;
  time_t currt = time(0);

  tmp = localtime(&currt);
  sprintf(timestring, "[%02d/%02d/%4d:%02d:%02d:%02d] ", 
	  tmp->tm_mday,tmp->tm_mon,1900 + tmp->tm_year, 
	  tmp->tm_hour,tmp->tm_min, tmp->tm_sec);
  return timestring;
}

void *yap_fullrealloc(void *ptr, size_t oldsize, size_t size) {
  void *p;
  if ( size<=oldsize ) {
    yap_message("yap_realloc() trying to shrink ... ignored\n");
    return ptr;
  }
  // yap_message("yap_realloc %x: %u - %u\n", ptr, oldsize, size);
  p = malloc(size);
  memset(p,0,size);
  memcpy(p,ptr,oldsize);
  free(ptr);
  return p;
}

/*
 * If logging to file, print time stamp, a message and optionally 
 * the system error string.  If syslogging, log a message and optionally 
 * the system error string.
 */
static void yap_print(int flag, int priority, const char *fmt, va_list ap)
{
  int errno_save;
  char buf[4096];

  assert(yap_fd>0 || use_syslog);
  buf[0] = 0;

  errno_save = errno;         /* value caller might want printed   */
  if ( !use_syslog && !(flag&YAP_MSG))
    strcpy(buf, yap_stamp());  /* prepend a message with time stamp */
  if (errno_save && (flag&YAP_ERRNO))
    sprintf(buf + strlen(buf), "%s: ", strerror(errno_save));
  vsprintf(buf + strlen(buf), fmt, ap);
  if ( use_syslog ) {
    syslog(priority, "%s", &buf[0]);
  } else {
    if ( write(yap_fd, buf, strlen(buf))<0 )
      exit(1);
  }
  errno = errno_save;
}


/*
 * Nonfatal error related to a system call.
 */
void yap_sysreport(const char *fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);
  yap_print(YAP_ERRNO, LOG_WARNING, fmt, ap);
  va_end(ap);
}

/*
 * Nonfatal error unrelated to a system call.
 */
void yap_report(const char *fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);
  yap_print(0, LOG_WARNING, fmt, ap);
  va_end(ap);
}

/*
 * Nonfatal message with no extra errno/stamp wanted
 */
void yap_message(const char *fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);
  yap_print(YAP_MSG, LOG_DEBUG, fmt, ap);
  va_end(ap);
}

void yap_va(const char *fmt, va_list ap)
{
  yap_print(YAP_MSG, LOG_DEBUG, fmt, ap);
}

/*
 * Fatal error related to a system call.
 */
void yap_sysquit(const char *fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);
  yap_print(YAP_ERRNO, LOG_CRIT, fmt, ap);
  va_end(ap);
  exit(1);
}

/*
 * Fatal error unrelated to a system call.
 */
void yap_quit(const char *fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);
  yap_print(0, LOG_CRIT, fmt, ap);
  va_end(ap);
  exit(1);
}

/*
 *  print the command line to the report file
 */
void yap_commandline(int argc, char**argv) {
  int loop;
  char buf[4096]; 
  sprintf(buf,"COMMAND-LINE: %s ", basename(argv[0])); 
  for ( loop=1; loop<argc; loop++ ) {
    sprintf(buf+strlen(buf),"%s ", argv[loop]);
  }
  strcat(buf,"\n");
  yap_report(buf);
}

char *yap_makename(char *stem, char *ext) {
  char *savefile;
  if ( !ext )
    ext = "";
  savefile = malloc(strlen(stem)+1+strlen(ext));
  if ( !savefile ) 
    yap_quit("Cannot allocate filename string\n");
  strcpy(savefile,stem);
  strcat(savefile,ext);
  return savefile;
}

int  yap_fileexists(char *fname) {
  struct stat buf;
  if ( stat(fname,&buf) )
    return 0;
  if ( S_ISREG(buf.st_mode) )
    return 1;
  return 0;
}
