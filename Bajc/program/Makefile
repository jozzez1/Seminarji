CXX=		clang++
PRGNAM=		test
CPUOPT= 	-march=corei7 -mssse3 -msse4.1
OFLAG=		-O3
OPTFLAGS=	${OFLAG} ${CPUOPT}
WARN=		-Wall
LIBFLAGS=	-L/usr/local/lib -lm
CXXFLAGS=	-I/usr/local/include/eigen3 \
		-I/usr/local/include \
		${WARN} ${OPTFLAG}

all:
	${CXX} ${CXXFLAGS} ${LIBFLAGS} main.cxx -o ${PRGNAM}
