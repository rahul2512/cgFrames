CC=g++
CFLAGS= -std=c++11 -fPIC -fpermissive
INCLUDES=-I./include -I$(NETCDF_ROOT)/include/ -I./netcdf-cxx4-4.3.0/include/ -I./armadillo-9.100.5/include/ 
#LIBS=-lnetcdf_c++4 -larmadillo 
LIBS=-L./netcdf-cxx4-4.3.0/lib -lnetcdf_c++4 -L./armadillo-9.100.5/ -larmadillo $(OPENBLAS_ROOT)/lib/libopenblas.a -lgfortran 
IDIR=include
SDIR=src
BDIR=bin
ODIR=build

_DEPS=atom_group.hpp PDBReader.hpp TopReader.hpp TrajReader.hpp AtomGroupReader.hpp Trajectory.hpp cgFrames.hpp utils.hpp fraWriter.hpp fitter_arma.hpp
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ=PDBReader.o TopReader.o TrajReader.o AtomGroupReader.o Trajectory.o cgFrames.o utils.o fraWriter.o fitter_arma.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


$(BDIR)/genFrames: $(SDIR)/genFrames.cpp $(OBJ)
	$(CC) $(CFLAGS) $(INCLUDES) $^ $(LIBS) -o $@

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CC) $(CFLAGS) $(INCLUDES) $< $(LIBS) -c -o $@

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o
