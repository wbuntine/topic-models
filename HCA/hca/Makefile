#  these should be inherited from above
# STB_HOME
# CFLAGS
# LDFLAGS 
# CPPFLAGS

CFILES = stats.c cache.c data.c gibbs.c like.c likesub.c lrs.c probs.c \
  samplea.c sampleb.c merge.c  \
  topics.c pctl.c sparsemap.c tprob.c phi.c \
   change.c  query.c loadhdp.c check.c 
HFILES = data.h hca.h probs.h sample.h stats.h pctl.h diag.h check.h \
	change.h

all:    hca

hca:	hca.o $(subst .c,.o,$(CFILES) ) ../util/libhca.a
	cc  $(CPPFLAGS) $(CFLAGS)  -o hca hca.o $(subst .c,.o,$(CFILES)) $(LDFLAGS)

../util/libhca.a:

clean: 
	rm -f *.o *.d hca gmon.out 

distclean:   clean
	rm -f hca[A-Z] hca[A-Z][A-Z] TAGS *.d

%.d: %.c
	@set -e; rm -f $@; \
	$(CC) -MM $(CPPFLAGS) $< > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

include $(CFILES:.c=.d)

