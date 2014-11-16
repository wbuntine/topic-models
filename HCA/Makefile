#
#  to get multi-threaded version working, add -DH_THREADS
#  to CPPFLAGS to Linux/Mac-OSX versions below
#

###########################################
#     installed location of where to find libstb
#       i.e.,   libstb.a is in $(STB)/lib/
STB_HOME := /home/wbuntine/
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
CPPFLAGS =  -I../util/ -I$(STB_HOME)/include/ $(MYDEBUG)
CFLAGS = -Wall 
LDFLAGS = -L../util -L$(STB_HOME)/lib -lhca -lstb$(MYLIB)  -lm -pthread

####################################################
##  Mac OSX versions
#CFLAGS = -D__WORDSIZE=64 -I/opt/local/include -Wall
#LDFLAGS = -L/opt/local/lib -L../util -L../lib -lhca -lstb -lm -pthread
#CPPFLAGS=-I/opt/local/include -I../util/ -I../lib/ 

export CFLAGS
export LDFLAGS
export CPPFLAGS

all:    util hca # tca

yes:

util:    yes
	cd util; make 

terms:    yes
	cd util; make terms

heap:    yes
	cd util; make heap

hca:    yes

hca:    yes
	cd hca; make 

tca:    yes
	cd tca; make 

clean:  yes
	cd util; make clean
	cd hca; make clean	
	cd tca; make clean	

distclean: 
	rm -f TAGS */TAGS
	cd util; make distclean
	cd hca; make distclean	
	cd tca; make distclean	

etags:  
	etags  util/*.c util/*.h dca/*.h cca/*.h
