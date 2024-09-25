# makefile for gaffer developed on Richard's Mac

#CFLAGS= -DMACOS -O3
CFLAGS= -DMACOS -g -target arm64-apple-macos11	# for debugging
#CFLAGS= -DLINUX -O3 -DOMP -fopenmp		# for OMP parallelisation - doesn't compile on Mac

ALL=gaffer seqconvert seqstat ONEview ONEstat syng syngprune seqextract

DESTDIR=~/bin

all: $(ALL)

install:
	cp $(ALL) $(DESTDIR)

clean:
	$(RM) *.o *~ $(ALL)
	\rm -r *.dSYM

### object files

UTILS_OBJS=hash.o dict.o array.o utils.o
UTILS_HEADERS=utils.h array.h dict.h hash.h
$(UTILS_OBJS): utils.h $(UTILS_HEADERS)

SEQIO_OPTS = -DONEIO
SEQIO_LIBS = -lm -lz 

seqio.o: seqio.c seqio.h 
	$(CC) $(CFLAGS) $(SEQIO_OPTS) -c $^

seqhash.o: seqhash.c seqhash.h 
	$(CC) $(CFLAGS) -c $^

kmerhash.o: kmerhash.c kmerhash.h 
	$(CC) $(CFLAGS) -c $^

ONElib.o: ONElib.c ONElib.h 
	$(CC) $(CFLAGS) -c $^

### programs

seqconvert: seqconvert.c seqio.o ONElib.o utils.o
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

seqextract: seqextract.c seqio.o ONElib.o dict.o array.o utils.o
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

seqstat: seqstat.c seqio.o ONElib.o array.o utils.o
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

gaffer: gaffer.c seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

syng: syng.c seqio.o seqhash.o kmerhash.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

syngprune: syngprune.c seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

ONEview: ONEview.c ONElib.o
	$(CC) $(CFLAGS) $^ -o $@ -lz

ONEstat: ONEstat.c ONElib.o
	$(CC) $(CFLAGS) $^ -o $@ -lz

### end of file
