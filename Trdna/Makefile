CC=gcc

SOURCE=integral.c
PRGNAM=calc

OFLAG=-O3
LIBFLAGS=-lm
CPUFLAGS=-march=native -mssse3 -msse4.2
OPTFLAGS=$(OFLAG) $(CPUFLAGS) -fstack-protector -funroll-loops -ffast-math

CFLAGS=-Wall $(OPTFLAGS) $(LIBFLAGS)

all:
	$(CC) $(SOURCE) -o $(PRGNAM) $(CFLAGS)

clean:
	rm -rf test calc output.txt
