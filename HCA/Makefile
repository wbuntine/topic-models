#     installed location of where to find libstb
#       i.e.,   libstb.a is in $(STB)/lib/
STB_HOME := /home/wbuntine/

##  Linux versions
#         add "-g" or "-O5 -DNDEBUG" or ....
#  Add "-DH_THREADS" to get multi-core version
#  WARNING:  threads and debugging code ('-g') don't mix because some of the
#            asserts are shuffled around by the compiler
# MYDEBUG=-O5 -DNDEBUG -ffast-math -DHAVE_INLINE -DH_THREADS
MYDEBUG=-g 
#  add "-g" for debugging libraries
MYLIB=
CPPFLAGS =  -I../util/ -I$(STB_HOME)/include/ $(MYDEBUG)
CFLAGS = -Wall 
LDFLAGS = -L../util -L$(STB_HOME)/lib -lddca -lstb$(MYLIB)  -lm -pthread

##  Mac OSX versions
#CFLAGS = -D__WORDSIZE=64 -I/opt/local/include -Wall -g 
#LDFLAGS = -L/opt/local/lib -lstb -lm
#CPPFLAGS=-I/opt/local/include 

export STB_HOME
export CFLAGS
export LDFLAGS
export CPPFLAGS

all:    util hca tca

yes:

util:    yes
	cd util; make 

terms:    yes
	cd util; make terms

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
