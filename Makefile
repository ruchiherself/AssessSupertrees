#uncomment first four lines for windows and comment next 8 lines
cpp=c++ -g -O -ansi -pedantic -Wall -W -Wno-long-long
cpp=c++ -g -O -funroll-loops -Wno-long-long
cc=gcc -O3 -fomit-frame-pointer -funroll-loops
OUTEXEC=AssessSupertree

# Mac OS X
# MAC_UNIVERSAL=-arch i386 -arch ppc -mmacosx-version-min=10.0
# cpp=c++ -g -O3 -fomit-frame-pointer -funroll-loops   

#${MAC_UNIVERSAL}
#cc=cc -O3 -fomit-frame-pointer -funroll-loops ${MAC_UNIVERSAL}
# OUTEXEC=AssessSupertree.macosx

# cpp=c++ -g -O3 -funroll-loops -Wno-long-long -pg
# cc=gcc -O3 -funroll-loops -pg

INCLUDE=-I./include
LIBRARY=

all: AssessSupertree
	
AssessSupertree: main.o rmq.o
	${cpp} main.o rmq.o ${INCLUDE} ${LIBRARY} -o ${OUTEXEC}

main.o: main.cpp Makefile 
	${cpp} ${INCLUDE} -c $<

rmq.o: rmq.c rmq.h Makefile
	${cc} -c $<

clean:
	rm -f *.o *~ core ${OUTEXEC}
