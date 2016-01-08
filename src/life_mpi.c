#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

int BS;
#define cell(_i_,_j_) board[ldboard*(_j_)+(_i_)]
#define ngb(_i_,_j_) nbngb[ldnbngb*((_j_)-1)+((_i_)-1)]

void output_board(int N, int *board, int ldboard, int loop)
{
    int i,j;
    printf("loop %d\n", loop);
    for (i=0; i<N; i++) {
	for (j=0; j<N; j++) {
	    if ( cell( i, j ) == 1)
		printf("X");
	    else
		printf(".");
	}
	printf("\n");
    }
}

int generate_initial_board(int N, int *board, int ldboard)
{
    int i, j, num_alive = 0;

    for (i = 0; i < N; i++) {
	for (j = 0; j < N; j++) {
	    if (i == N/2 || j == N/2) {
		cell(i, j) = 1;
		num_alive ++;
	    }
	    else {
		cell(i, j) = 0;
	    }
	}
    }

    return num_alive;
}

// Synchronized implementation
int main(int argc, char *argv[]){
  int i, j, loop, num_alive, maxloop;
  int ldboard, ldnbngb;
  double t1, t2;
  double temps;
  int size, rank; 
  
  
  int *board;
  int *nbngb;
  
  if (argc < 3) {
    printf("Usage: %s nb_iterations size\n", argv[0]);
    return EXIT_SUCCESS;
  } else {
    maxloop = atoi(argv[1]);
    BS = atoi(argv[2]);
    //printf("Running sequential version, grid of size %d, %d iterations\n", BS, maxloop);
  }
  num_alive = 0;

  /* Leading dimension of the board array */
  ldboard = BS + 2;
  /* Leading dimension of the neigbour counters array */
  ldnbngb = BS;

  board = malloc( ldboard * ldboard * sizeof(int) );
  nbngb = malloc( ldnbngb * ldnbngb * sizeof(int) );

  num_alive = generate_initial_board( BS, &(cell(1, 1)), ldboard );

  printf("Starting number of living cells = %d\n", num_alive);
  
  MPI_Init(NULL, NULL);
  MPI_Comm_size(&size, MPI_COMM_WORLD);
  MPI_Comm_rank(&rank, MPI_COMM_WORLD);
  
  
  return EXIT_SUCCESS;
}