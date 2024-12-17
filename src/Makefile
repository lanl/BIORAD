OBJS = options.o mpi_util.o parse_util.o pdb.o cluster.o regression.o features.o seq_overlap.o \
	rosetta.o hbond.o linear_regression.o

CC = mpic++
cc = mpicc

FLAGS = -O3 -Wall -Wno-multichar

PROFILE = #-g #-pg
OPENMP = -fopenmp
	
INC = -I. -I$(HOME)/gsl/include
LIBS = -lm $(HOME)/gsl/lib/libgsl.a $(HOME)/gsl/lib/libgslcblas.a

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


