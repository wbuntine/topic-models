#  these should be inherited from above
# STB_HOME
# CFLAGS
# LDFLAGS 
# CPPFLAGS

CFILES = stats.c cache.c data.c gibbs.c like.c likesub.c lrs.c probs.c \
  samplea.c sampleb.c topics.c pctl.c  change.c
HFILES = data.h tca.h probs.h sample.h stats.h pctl.h change.h

all:    tca

tca:	tca.o $(subst .c,.o,$(CFILES) ) ../util/libutil.a
	cc  $(CPPFLAGS) $(CFLAGS)  -o tca tca.o $(subst .c,.o,$(CFILES)) $(LDFLAGS)

../util/libutil.a:

clean: 
	rm -f *.o *.d tca gmon.out 

distclean:   clean
	rm -f tca[A-Z] tca[A-Z][A-Z] TAGS *.d

%.d: %.c
	@set -e; rm -f $@; \
	$(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

include $(CFILES:.c=.d)

