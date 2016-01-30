CC = mpicc
WFLAGS = -Wall -Wextra
MATHS = -lm
OMP = -fopenmp
EXE = exe_omp exe_seq exe_mpi
MPIR = mpirun -np 4

.PHONY: all

all: $(EXE)

exe_mpi: life_mpi.o
	$(CC) -o $@ $^ $(WFLAGS) $(MATHS)

exe_omp: life_openmp.o
	 $(CC) -o $@ $^ $(WFLAGS) $(OMP)

exe_seq: life_seq.o
	$(CC) -o $@ $^ $(WFLAGS) $(OMP)

%.o: %.c
	$(CC) -c $^ $(WFLAGS) $(OMP)
clean:
	rm *.o
	rm $(EXE)
run:
	./exe_seq
	./exe_omp
	$(MPIR) -h
