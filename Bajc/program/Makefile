CC=clang
PRGNAM=fermion
SOURCE=fitter.c

.if !empty (OSTYPE:Mlinux*) || defined (FORCE_CC)
  CC=cc
.endif

OFLAG=-O3
CPUFLAGS=-march=corei7 -mssse3 -msse4.1 
OPTFLAGS=$(OFLAG) -funroll-loops -fno-strict-aliasing -ffast-math -fstack-protector
WARNFLAGS=-Wall -pedantic-errors
LIBFLAGS=-lgsl -lgslcblas -lm -lc -I/usr/local/include -L/usr/local/lib

CFLAGS=$(CPUFLAGS) $(WARNFLAGS) $(LIBFLAGS) $(OPTFLAGS)

all: runner
	$(CC) $(CFLAGS) -g -o $(PRGNAM) $(SOURCE) 

runner:
	$(CC) $(CFLAGS) -g -o drsi runner.c
