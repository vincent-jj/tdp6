CC = gcc
WFLAGS = -Wall -Wextra
OMP = -fopenmp
EXE = exe_omp exe_seq

.PHONY: all

all: $(EXE)

exe_omp: life_openmp.o
	 $(CC) -o $@ $^ $(WFLAGS) $(OMP)

exe_seq: life_seq.o
	$(CC) -o $@ $^ $(WFLAGS) $(OMP)

%.o: %.c
	$(CC) -c $^ $(WFLAGS) $(OMP)
clean:
	rm *.o
	rm $(EXE)
