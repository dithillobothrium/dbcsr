include ../../arch/local.mdbg


OBJ        = carma_mm.o pdgemm_dbcsr.o 

EXE        = pdgemm_dbcsr
INCLUDES   += -I../../obj/local/mdbg
LIBS      :=  -lcp2kmpiwrap -lcp2kbase \
	     -ldbcsr -ldbcsrmm -ldbcsrops -ldbcsrbblock -ldbcsrbdist -ldbcsrdata \
             -lcp2kacc -ldbcsrbase  -lstdc++ -L../../lib/local/mdbg $(LIBS)
#LIBS      +=
#LIBS      += /project/s238/alazzaro/libsmm/libxsmm/sandybridge/lib/libsmm_dnn_cray.gnu.a

FLAGS = -O0 -g -fopenmp -ffree-form -fimplicit-none -std=f2003

%.o: %.F
	$(FC) $(FLAGS) $(INCLUDES)  $< $(LIBS) -c -o $@

$(EXE): $(OBJ)
	$(FC)  $(FLAGS) $(INCLUDES)  $^ $(LIBS) -o  $@

clean:
	rm -f *~ $(EXE) *.o *.mod
