OBJS = options.o mpi_util.o parse_util.o pdb.o cluster.o regression.o features.o seq_overlap.o \
	rosetta.o hbond.o linear_regression.o

FLAGS = -O3 -Wall -Wno-multichar -std=c++17

# Uncomment to enable parallel processing with MPI
#USE_MPI = 1
ifdef USE_MPI
	CC = mpic++
	cc = mpicc
	FLAGS += -DUSE_MPI
else
	CC = c++
	cc = cc
endif

PROFILE = #-g #-pg
OPENMP = #-fopenmp
	
INC = -I. 
LIBS = -lm

# Uncomment to enable the calculation of linear regressions with the Gnu Scientific Library
#USE_GSL = 1
ifdef USE_GSL
	FLAGS += -DUSE_GSL
	INC += -I$(HOME)/gsl/include
	LIBS += $(HOME)/gsl/lib/libgsl.a $(HOME)/gsl/lib/libgslcblas.a
endif

.SUFFIXES : .o .cpp .c
.cpp.o:
	$(CC) $(FLAGS) $(PROFILE) $(OPENMP) $(INC) -c $<
 
.c.o:
	$(cc) $(FLAGS) $(PROFILE) $(OPENMP) $(INC) -c $<

all: biorad
	 	
biorad: $(OBJS) main.o 
	$(CC) $(PROFILE) -o biorad $(OBJS) $(OPENMP) main.o $(LIBS)

clean:
	-rm -f *.o


