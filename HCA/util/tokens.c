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
 *    read tokens in range W0:WE-1
 *    place in vector offset by W0
 *    allocate all string memory
 */
#define MAXSTRING 1000
char **read_vocab(char *infile, int W0, int WE, int len) {
  /*
   *   set up token names in vector
   */
  int W = WE-W0;
  char *vocfile = malloc(strlen(infile)+10);
  int i;
  char *cmem = (char *)malloc(len*W+MAXSTRING);
  int memused = 0;
  int memallocd = len*W;
  int *cvec = (int *)malloc(W*sizeof(int));
  char **ctmp = (char **)malloc(W*sizeof(char *));
  FILE *fr;
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
  for (i=0;i<W0;i++) {
    if ( fgets(&cmem[cvec[0]],MAXSTRING-1,fr)==NULL )
      yap_sysquit("Cannot read line %d/%d from '%s'\n", i+1, W, vocfile);
  }
  for (i=0;i<W;i++) {
    int sl;
    if ( fgets(&cmem[cvec[i]],MAXSTRING-1,fr)==NULL )
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
      if ( memused+MAXSTRING>=memallocd ) {
	int reall=(W-i)*len;
	assert(memused<memallocd);
	if ( memused+MAXSTRING>=reall+memallocd )
	  reall = MAXSTRING;
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
