
EXE        = pdgemm_dbcsr
INCLUDES   = -I/project/s238/isivkov/dbcsr/obj/CRAY-XC50-gfortran-cuda/psmp
LIBS       = -L/project/s238/isivkov/dbcsr/lib/CRAY-XC50-gfortran-cuda/psmp -lcp2kmpiwrap -lcp2kbase \
	     -ldbcsr -ldbcsrmm -ldbcsrops -ldbcsrbblock -ldbcsrbdist -ldbcsrdata \
             -lcp2kacc -ldbcsrbase -lcp2kacccuda -lcusmm -lstdc++
#LIBS      +=
LIBS      += /project/s238/alazzaro/libsmm/libxsmm/sandybridge/lib/libsmm_dnn_cray.gnu.a

$(EXE): $(EXE).F
	ftn -O0 -g $(INCLUDES) -fopenmp -ffree-form -fimplicit-none -std=f2003 $< $(LIBS) -o $@

clean:
	rm -f *~ $(EXE)
