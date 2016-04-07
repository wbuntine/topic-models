#
#  to get multi-threaded version working, add -DH_THREADS
#  to CPPFLAGS to Linux/Mac-OSX versions below
#

###########################################
#     Alter these to setup your local
#     directory structure, ...
#     but the defaults should work

#     top level for HCA
PWD := $(shell pwd)
#     installed location of where to find libstb
#       i.e.,   libstb.a is in $(STB)/lib/
STB_HOME := $(PWD)/../../libstb
export STB_HOME
export PWD

#######################################
##  Linux versions
#
#  add "-g" or "-O5 -DNDEBUG" or ....
#  WARNING:  threads and debugging code ('-g') don't mix because some of the
#            asserts are shuffled around by the compiler
# MYDEBUG=-pg -DNDEBUG -ffast-math -DHAVE_INLINE
# MYDEBUG=-O5 -DNDEBUG -ffast-math -DHAVE_INLINE # -DH_THREADS
MYDEBUG=-g
CPPFLAGS :=  -I$(PWD)/util/ -I$(STB_HOME)/lib/ $(MYDEBUG)
CFLAGS = -Wall
LDFLAGS := -L$(PWD)/util -L$(STB_HOME)/lib
LIBS = -lhca -lstb -lm -pthread

CC = gcc

####################################################
##  Mac OSX versions
#CFLAGS = -D__WORDSIZE=64 -I/opt/local/include -Wall
#LDFLAGS = -L/opt/local/lib -L../util -L../lib -lhca -lstb -lm -pthread
#CPPFLAGS=-I/opt/local/include -I../util/ -I../lib/

export CC
export CFLAGS
export LDFLAGS
export LIBS
export CPPFLAGS

all:    util hca tca

yes:

util:    yes
	cd util; ${MAKE}

terms:    yes
	cd util; ${MAKE} terms

heap:    yes
	cd util; ${MAKE} heap

hca:    yes

hca:    yes
	cd hca; ${MAKE}

tca:    yes
	cd tca; ${MAKE}

clean:  yes
	cd util; ${MAKE} clean
	cd hca; ${MAKE} clean
	cd tca; ${MAKE} clean

distclean:
	rm -f TAGS */TAGS
	cd util; ${MAKE} distclean
	cd hca; ${MAKE} distclean
	cd tca; ${MAKE} distclean

etags:
	etags  util/*.c util/*.h dca/*.h cca/*.h
