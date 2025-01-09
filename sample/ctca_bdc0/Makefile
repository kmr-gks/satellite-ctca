PREFIX=./ctca
#PREFIX=/home/tmp/cotocoa/cotocoa-1.2.3-intel2023.2
CC=mpiicc
FC=mpiifort
FFLAGS=-I${PREFIX}/include
LDFLAGS=-L${PREFIX}/lib
LIBS=-lctca
F_LIBS=-lctca_f -lctca_mod

PROGS = sample-bdc0_requester sample-bdc0_coupler sample-bdc0_worker1

all: ${PROGS}

sample-bdc0_requester: sample-bdc0_requester.f90
	$(FC) $(FFLAGS) $(LDFLAGS) sample-bdc0_requester.f90 $(F_LIBS) $(LIBS) -o sample-bdc0_requester

sample-bdc0_coupler: sample-bdc0_coupler.f90
	$(FC) $(FFLAGS) $(LDFLAGS) sample-bdc0_coupler.f90 $(F_LIBS) $(LIBS) -o sample-bdc0_coupler

sample-bdc0_worker1: sample-bdc0_worker1.f90
	$(FC) $(FFLAGS) $(LDFLAGS) sample-bdc0_worker1.f90 $(F_LIBS) $(LIBS) -o sample-bdc0_worker1

clean:
	/bin/rm -f *.o *.mod ${PROGS}
