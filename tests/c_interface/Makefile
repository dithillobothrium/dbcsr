
EXE = dbcsr_example_1_c dbcsr_example_1_cpp

DBCSR      = /data/isivkov/dbcsr_git
INCLUDES  += -I$(DBCSR)/obj/local/mdbg
LIBS      += -L$(DBCSR)/lib/local/mdbg \
	     -ldbcsr -ldbcsrmm -ldbcsrops -ldbcsrbblock -ldbcsrbdist -ldbcsrdata \
	     -lcp2kacc -ldbcsrbase -lcp2kmpiwrap -lcp2kbase
LIBS      += -L${SCALAPACK_PATH}/lib -lscalapack -lreflapack -lrefblas

MPI_LINK_FLAGS = $(shell mpif90 --showme:link)

default: $(EXE)

#dbcsr_example_1_c: dbcsr_example_1.c dbcsr.h
#	mpicc -O3 -c $<
#	mpicc -O3 -fopenmp $(<:.c=.o) -lmpi_mpifh -lgfortran -lm  $(LIBS) -o $@

dbcsr_example_1_cpp: dbcsr_example_1.cpp dbcsr.h
	mpicxx -O3 -fopenmp -std=c++11 -c $<
	mpicxx -O3 -fopenmp -std=c++11 $(<:.c=.o) $(LIBS) -lmpi_mpifh  -lgfortran -o $@

clean:
	rm -rf *~ $(EXE) *.o
