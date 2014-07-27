#
#   CFLAGS and LDFLAGS should be in environment vars

SRC = beta.c coher.c csparse.c myarms.c tokens.c util.c yap.c pmi.c \
	dread.c misi.c  ehash.c fvec.c getline.c heap.c
HFILES = myarms.h util.h yap.h dread.h misi.h pmi.h ehash.h fvec.h heap.h
OBJ = $(SRC:.c=.o)
LIBRARY = libhca.a

library:  $(OBJ)
	ar rcs $(LIBRARY) $(OBJ)

terms:  terms.c
	$(CC) -g -o terms terms.c util.o yap.o dread.o tokens.o -lm

heap:  heap.c
	$(CC) -DHEAP_MAIN -g -o heap heap.c -lm

clean: 
	rm -f *.o $(LIBRARY) terms heap

distclean:

etags:  $(SRC) $(HFILES)
	etags  $(SRC) $(HFILES)

depend:
	makedepend $(FILES) 

# DO NOT DELETE THIS LINE -- make depend depends on it. 
