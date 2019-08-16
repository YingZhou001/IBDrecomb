IDIR=.
CC=gcc
CFLAGS=-I$(IDIR)

LIBS=-lm

OBJ=IBDrecomb.o check.o read.o tools.o end.o estimateRec.o

DEPS=head.h check.h read.h tools.h end.h estimateRec.h

IBDrecomb:$(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

IBDrecomb.o:IBDrecomb.c $(DEPS)
	${CC} -c -o $@ $< $(CFLAGS)

check.o:check.c $(DEPS)
	${CC} -c -o $@ $< $(CFLAGS)

end.o:end.c $(DEPS)
	${CC} -c -o $@ $< $(CFLAGS)

read.o:read.c $(DEPS)
	${CC} -c -o $@ $< $(CFLAGS)

estimateRec.o:estimateRec.c $(DEPS)
	${CC} -c -o $@ $< $(CFLAGS)

tools.o:tools.c $(DEPS)
	${CC} -c -o $@ $< $(CFLAGS)


.PHONY:clean

clean:
	rm -f *.o
