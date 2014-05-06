/*
 * Load vocabulary
 * Copyright (C) 2011 Wray Buntine 
 * All rights reserved.
 *
 * This Source Code Form is subject to the terms of the Mozilla 
 * Public License, v. 2.0. If a copy of the MPL was not
 * distributed with this file, You can obtain one at
 *      http://mozilla.org/MPL/2.0/.
 *
 * Author: Wray Buntine (wray.buntine@nicta.com.au)
 *
 *     input filename and read vocab;
 *     assumes all names < 30 chars!!
 *
 *  result can be freed with a single free()
 *
 *  horrible hacks, so use for testing only
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include "assert.h"
#include "yap.h"

/*
 *    len is just a hint, it will resize if more needed;
 *    to free result:
 *          free(ctmp[0]); free(ctmp);
 */
char **read_vocab(char *infile, int W, int len) {
  /*
   *   set up token names in vector
   */
  char *vocfile = malloc(strlen(infile)+10);
  int i;
  char *cmem = (char *)malloc(len*W);
  int memused = 0;
  int memallocd = len*W;
  int *cvec = (int *)malloc(W*sizeof(int));
  char **ctmp = (char **)malloc(W*sizeof(char *));
  FILE *fr;
  char wordbuf[1000];
  if ( !ctmp || !cmem ) 
    yap_quit("Cannot allocate name space in read_vocab()\n");
  /*
   *  open file
   */
  strcpy(vocfile,infile);
  if ( strrchr(vocfile,'.') )
    strcpy(strrchr(vocfile,'.'),".tokens");
  else
    strcpy(vocfile+strlen(vocfile),".tokens");
  fr = fopen(vocfile,"r");
  if ( !fr ) 
    yap_sysquit("vocfile '%s' not opened\n", vocfile);

  /*
   *  read terms
   */  
  cvec[0] = 0;
  for (i=0;i<W;i++) {
    int sl;
    if ( fgets(&cmem[cvec[i]],sizeof(wordbuf)-1,fr)==NULL )
      yap_sysquit("Cannot read line %d/%d from '%s'\n", i+1, W, vocfile);
    sl = strlen(&cmem[cvec[i]]);
    if ( ! iscntrl(cmem[cvec[i]+sl-1]) ) {
      /*
       *   line too long
       */
      yap_quit("Cannot parse line %d/%d from '%s', too long\n", i+1, W, vocfile);
    }
    /*
     *  assumes whole line read
     */
    while ( sl>0 && iscntrl(cmem[cvec[i]+sl-1]) ) {
      cmem[cvec[i]+sl-1] = 0;
      sl--;
    }
    memused += sl+1;
    if ( i<W-1 ) {
      cvec[i+1] = cvec[i]+sl+1;
      if ( memused+sizeof(wordbuf)>=memallocd ) {
	int reall=(W-i)*len;
	assert(memused<memallocd);
	if ( memused+sizeof(wordbuf)>=reall+memallocd )
	  reall = sizeof(wordbuf);
	memallocd += reall;
	cmem = realloc(cmem, memallocd);
	if ( !cmem )     
	  yap_quit("Cannot reallocate name space in read_vocab()\n");

      }
    }
  }
  fclose(fr);

  /*
   *   now form pointers;
   *   could not do earlier since mem may be reallocd
   */
  for (i=0; i<W; i++)
    ctmp[i] = cmem + cvec[i];
  free(cvec);

  free(vocfile);
  return ctmp;
}

void free_vocab(char **ctmp) {
	 free(ctmp[0]); free(ctmp);
}
