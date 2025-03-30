#
# Makefile for snarkhunter
#

SHELL = /bin/sh
CC=gcc

# Compiling executing program with DWORDSIZE=32 is slightly faster, 
# but limits the order of the graphs to 32.
# -mpopcnt is necessary for __builtin_popcount
CC32 = ${CC} -DWORDSIZE=32 -DMAXN=WORDSIZE -march=native
CC64 = ${CC} -DWORDSIZE=64 -DMAXN=WORDSIZE -march=native
CFLAGS = -O3

all : 
	rm -rf polehunter 
	${CC32} ${CFLAGS} polehunter.c nautyW1.a -o polehunter

64bit :
	rm -rf polehunter-64	
	${CC64} ${CFLAGS} polehunter.c nautyL1.a -o polehunter-64

profile :
	rm -rf polehunter-profile
	${CC32} -pg -g polehunter.c nautyW1.a -o polehunter-profile

edgepent :
	rm -rf edgepent
	${CC32} ${CFLAGS} -DPLUGIN='"edgepent.c"' polehunter.c nautyW1.a -o edgepent
