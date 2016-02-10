#
#  to get multi-threaded version working, add -DH_THREADS
#  to CPPFLAGS to Linux/Mac-OSX versions below
#

###########################################
#     installed location of where to find libstb
#       i.e.,   libstb.a is in $(STB)/lib/
STB_HOME := $(shell pwd)/../../libstb
export STB_HOME

#######################################
##  Linux versions
#
#  add "-g" or "-O5 -DNDEBUG" or ....
#  WARNING:  threads and debugging code ('-g') don't mix because some of the
#            asserts are shuffled around by the compiler
# MYDEBUG=-pg -DNDEBUG -ffast-math -DHAVE_INLINE # -DH_THREADS
# MYDEBUG=-O5 -DNDEBUG -ffast-math -DHAVE_INLINE -DH_THREADS
MYDEBUG=-g
CPPFLAGS :=  -I$(shell pwd)/util/ -I$(STB_HOME)/lib/ $(MYDEBUG)
CFLAGS = -Wall
LDFLAGS := -L$(shell pwd)/util -L$(STB_HOME)/lib
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
