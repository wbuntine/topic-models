#
#   CFLAGS and LDFLAGS should be in environment vars

SRC = beta.c coher.c csparse.c myarms.c tokens.c util.c yap.c pmi.c dread.c misi.c 
HFILES = myarms.h util.h yap.h dread.h misi.h
OBJ = $(SRC:.c=.o)
LIBRARY = libddca.a

library:  $(OBJ)
	ar rcs $(LIBRARY) $(OBJ)

terms:  terms.c
	$(CC) -g -o terms terms.c util.o yap.o dread.o tokens.o -lm

clean: 
	rm -f *.o $(LIBRARY) terms

distclean:

etags:  $(SRC) $(HFILES)
	etags  $(SRC) $(HFILES)

depend:
	makedepend $(FILES) 

# DO NOT DELETE THIS LINE -- make depend depends on it. 
