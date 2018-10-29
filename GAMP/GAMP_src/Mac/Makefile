# makefile for gamp

SRC     = .

OPTS    = -DTRACE -DENAGLO -DENAGAL -DENACMP -DENAQZS -DNFREQ=3
#OPTS    = -DENAGLO -DENAQZS -DENAGAL -DENACMP -DNFREQ=2

# for no lapack
CC      = gcc
#CFLAGS  = -g -Wall -O3 -ansi -pedantic -Wno-unused-but-set-variable -I$(SRC) $(OPTS)
#CFLAGS  = -g -Wall -O3 -DDEBUG -std=gnu89 -pedantic -I$(SRC) $(OPTS)
CFLAGS  = -g -Wall -O3 -std=gnu89 -pedantic -I$(SRC) $(OPTS)
LDLIBS  = -lm -lrt

#CFLAGS  = -Wall -O3 -ansi -pedantic -Wno-unused-but-set-variable -I$(SRC) -DLAPACK $(OPTS)
#LDLIBS  = -lm -lrt -llapack -lblas

# for gprof
#CFLAGS  = -Wall -O3 -ansi -pedantic -Wno-unused-but-set-variable -I$(SRC) -DLAPACK $(OPTS) -pg
#LDLIBS  = -lm -lrt -llapack -lblas -pg

# for mkl
##MKLDIR  = /opt/intel/mkl
#MKLDIR  = /proj/madoca/lib/mkl
#CFLAGS  = -O3 -ansi -pedantic -Wno-unused-but-set-variable -I$(SRC) $(OPTS) -DMKL
#LDLIBS  = -L$(MKLDIR)/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lgomp -lm -lrt

all    : gamp
gamp   : gamp.o rtkcmn.o rinex.o gampPos.o results.o tides.o
gamp   : lambda.o preceph.o spp.o brdceph.o ppp.o ionex.o
gamp   : prePpp.o myPpp.o myStr.o myMath.o myRinex.o myRtkcmn.o

gamp.o     : $(SRC)/gamp.c
	$(CC) -c $(CFLAGS) $(SRC)/gamp.c
rtkcmn.o   : $(SRC)/rtkcmn.c
	$(CC) -c $(CFLAGS) $(SRC)/rtkcmn.c
rinex.o    : $(SRC)/rinex.c
	$(CC) -c $(CFLAGS) $(SRC)/rinex.c
gampPos.o  : $(SRC)/gampPos.c
	$(CC) -c $(CFLAGS) $(SRC)/gampPos.c
results.o  : $(SRC)/results.c
	$(CC) -c $(CFLAGS) $(SRC)/results.c
lambda.o   : $(SRC)/lambda.c
	$(CC) -c $(CFLAGS) $(SRC)/lambda.c
preceph.o  : $(SRC)/preceph.c
	$(CC) -c $(CFLAGS) $(SRC)/preceph.c
spp.o      : $(SRC)/spp.c
	$(CC) -c $(CFLAGS) $(SRC)/spp.c
brdceph.o  : $(SRC)/brdceph.c
	$(CC) -c $(CFLAGS) $(SRC)/brdceph.c
ppp.o      : $(SRC)/ppp.c
	$(CC) -c $(CFLAGS) $(SRC)/ppp.c
ionex.o    : $(SRC)/ionex.c
	$(CC) -c $(CFLAGS) $(SRC)/ionex.c
tides.o    : $(SRC)/tides.c
	$(CC) -c $(CFLAGS) $(SRC)/tides.c
prePpp.o   : $(SRC)/prePpp.c
	$(CC) -c $(CFLAGS) $(SRC)/prePpp.c
myPpp.o    : $(SRC)/myPpp.c
	$(CC) -c $(CFLAGS) $(SRC)/myPpp.c
myStr.o    : $(SRC)/myStr.c
	$(CC) -c $(CFLAGS) $(SRC)/myStr.c
myMath.o    : $(SRC)/myMath.c
	$(CC) -c $(CFLAGS) $(SRC)/myMath.c
myRinex.o    : $(SRC)/myRinex.c
	$(CC) -c $(CFLAGS) $(SRC)/myRinex.c
myRtkcmn.o    : $(SRC)/myRtkcmn.c
	$(CC) -c $(CFLAGS) $(SRC)/myRtkcmn.c

clean :
	rm -f *.o

