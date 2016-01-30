#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include <mpi.h>
#include <math.h>

#define cell( _i_, _j_ ) board[ ldboard * (_j_) + (_i_) ]
#define ngb( _i_, _j_ )  nbngb[ ldnbngb * ((_j_) - 1) + ((_i_) - 1 ) ]

int BS;

// Cell neighbour's composition
// 4 2 5
// 0 X 1
// 6 3 7

//Neighbors 
int L=0,R=1,U=2,D=3,UL=4,UR=5,LR=6,LL=7;

#define ROW 0
#define COLUMN 1

inline void helper(){
  printf("Program Option :\n\t* First arg : nb iterations in loop\n\t* Second arg : total number of cells\n\t* Help : -h");
  return;
}

inline double mytimer(void){
  struct timeval tp;
  gettimeofday( &tp, NULL );
  return tp.tv_sec + 1e-6 * tp.tv_usec;
}

void output_board(int N, int *board, int ldboard, int loop){
  int i,j;
  fprintf(stderr, "loop %d\n", loop);
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      if ( cell( i, j ) == 1)
	fprintf(stderr, "X");
      else
	fprintf(stderr, " ");
    }
    fprintf(stderr,"\n");
  }
}

/**
 * This function generates the iniatl board with one row and one
 * column of living cells in the middle of the board
 */
int generate_initial_board(int N, int *board, int ldboard){
  int i, j, num_alive = 0;
  for (i = 1; i <= N; i++) {
    for (j = 1; j <= N; j++) {
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

/*
 * This function process a nbngb block --> process frontier blocks
 */
void process_frontier(int k, int blocksize, int* board, int frontier, int ldboard, int* nbngb, int ldnbngb){
  /* Different process if shared from a column or block frontier */
  if (frontier == ROW){
    int i = k;
    int j;
    for (j = 1; j <= blocksize; ++j) {
      ngb( i, j ) =
	cell( i-1, j-1 ) + cell( i, j-1 ) + cell( i+1, j-1 ) +
	cell( i-1, j   ) +                  cell( i+1, j   ) +
	cell( i-1, j+1 ) + cell( i, j+1 ) + cell( i+1, j+1 );
    }
    return;
  }
  if (frontier==COLUMN){
    int j = k;
    int l;
    for (l = 1; l <= blocksize; ++l) {
      ngb( l, j ) =
	cell( l-1, j-1 ) + cell( l, j-1 ) + cell( l+1, j-1 ) +
	cell( l-1, j   ) +                  cell( l+1, j   ) +
	cell( l-1, j+1 ) + cell( l, j+1 ) + cell( l+1, j+1 );
    }
    return;
  }
  printf("Error, not compatible frontier type in process frontier.\nExiting program\n");
  MPI_Finalize();
  exit(EXIT_FAILURE);
}

void inter_proc_communications(MPI_Request* request, int* neighbours, MPI_Comm grid, int blocksize, int* board, int ldboard, MPI_Datatype block_line){
  // Sending / Receving informations to / from the left / right process
  MPI_Isend(&cell(1, 1), blocksize, MPI_INT, neighbours[L], 0, grid, &request[0]);
  MPI_Irecv(&cell(1, blocksize+1), blocksize, MPI_INT, neighbours[R], 0, grid, &request[0]);

  // Sending / Receving informations to / from the right / left process
  MPI_Isend(&cell(1, blocksize), blocksize, MPI_INT, neighbours[R], 0, grid, &request[1]);
  MPI_Irecv(&cell(1, 0), blocksize, MPI_INT, neighbours[L], 0, grid, &request[1]);     

  // Sending / Receving informations to / from the up-right / low-left process
  MPI_Isend(&cell(1, blocksize), 1, MPI_INT, neighbours[UR], 0, grid, &request[2]);
  MPI_Irecv(&cell(blocksize+1, 0), 1, MPI_INT, neighbours[LL], 0, grid, &request[2]); 

  // Sending / Receving informations to / from the low-right / up-left process
  MPI_Isend(&cell(blocksize, blocksize), 1, MPI_INT, neighbours[LR], 0, grid, &request[3]);
  MPI_Irecv(&cell(0, 0), 1, MPI_INT, neighbours[UL], 0, grid, &request[3]);

  // Sending / Receving informations to / from the up-left / low-right process
  MPI_Isend(&cell(1, 1), 1, MPI_INT, neighbours[UL], 0, grid, &request[4]);
  MPI_Irecv(&cell(blocksize+1, blocksize+1), 1, MPI_INT, neighbours[LR], 0, grid, &request[4]);

  // Sending / Receving informations to / from the low-left / up-right process
  MPI_Isend(&cell(blocksize, 1), 1, MPI_INT,neighbours[LL], 0, grid, &request[5]);
  MPI_Irecv(&cell(0, blocksize+1), 1, MPI_INT, neighbours[UR], 0, grid, &request[5]);

  // Sending / Receving informations to / from the bottom / up process
  MPI_Isend(&cell(blocksize, 1), 1, block_line,neighbours[D], 0, grid, &request[6]);
  MPI_Irecv(&cell(0, 1), 1, block_line, neighbours[U], 0, grid, &request[6]);

  // Sending / Receving informations to / from the up / bottom process
  MPI_Isend(&cell(1, 1), 1, block_line,neighbours[U], 0, grid, &request[7]);
  MPI_Irecv(&cell(blocksize+1, 1), 1, block_line, neighbours[D], 0, grid, &request[7]);
}

/* Function that creates a list of the neighbours for each processes */ 
void neighbour_table(int* neighbours, MPI_Comm grid, int proc_rank){
  int move, id;
  int coord[2];
  id = 1;
  // move from a proc to an other ==> move = 1
  move = 1;
  // get Left and Right
  MPI_Cart_shift(grid, id, move, &neighbours[L], &neighbours[R]);
  id = 0;
  // get Up and Down
  MPI_Cart_shift(grid, id, move, &neighbours[U], &neighbours[D]);
  // get current proc coordinates
  MPI_Cart_coords(grid, proc_rank, 2, coord);
  coord[0]--;
  coord[1]--;
  // determine Up-Left neighbour
  MPI_Cart_rank(grid, coord, &neighbours[UL]);
  coord[1]+=2;
  // determine Up-Right neighbour
  MPI_Cart_rank(grid, coord, &neighbours[UR]);
  coord[0]+=2;
  // determine Down-Right neighbour
  MPI_Cart_rank(grid, coord, &neighbours[LR]);
  coord[1]-=2;
  // determine Down-Left neighbour
  MPI_Cart_rank(grid, coord, &neighbours[LL]);
  return;
}

int main(int argc, char* argv[]){
  MPI_Init(NULL, NULL);
  int rank, size;
  int loop, num_alive, loop_iterations;
  int ldboard, ldnbngb, ldglobalboard;
  double t1, time, final_time;
  int periods[2] = {1, 1};
  int *globboard= NULL;
  int *globboard2= NULL;
  int *board;
  int *nbngb;

  /* Initialization of MPI */
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size);
  if(argc >= 2){
    if(!strcmp("-h",argv[1])){
      if(!rank)
	helper();
      MPI_Finalize();
      return EXIT_SUCCESS;
    }
  }
  int i, j;
  int process_per_row = sqrt(size);
  int process_per_column = sqrt(size);
  int dims[2] = {process_per_row, process_per_column};
  
  // It only works if the number of process in the input is a perfect square
  if(size != process_per_column*process_per_row){
    fprintf(stderr, "Square Perfect needed as input size.\nExiting Program.");
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  MPI_Comm grid;

  // Initialize cartesian grid
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods,0, &grid);
  MPI_Comm_rank(grid, &rank);

  /* User input */
  if (argc < 2) {
    loop_iterations = 10;
    BS = 30;
  } else if (argc >= 2){
    loop_iterations = atoi(argv[1]);
    if(argc > 2)
      BS = atoi(argv[2]);
    else
      BS = 30;
  }
  num_alive = 0;

  /*Leading dimension of global board array*/
  ldglobalboard = BS + 2; // +2 because of upper and above added (+ X +)
  /* Leading dimension of board array */
  ldboard = BS/process_per_row + 2; // +2 because of upper and above added (+ X +)
  /* Leading dimension of neigbour array */
  ldnbngb = BS/sqrt(size); // Same number of element in each process which is equal to this formula

  // Initialization of cells board
  board = (int *)malloc( ldboard * ldboard * sizeof(int) );
  nbngb = (int *)malloc( ldnbngb * ldnbngb * sizeof(int) );

  // Initialization of global cell board (which is common between all processes)
  if(!rank){
    globboard = (int *)malloc(ldglobalboard*ldglobalboard * sizeof(int));
    globboard2 = (int *)malloc(ldglobalboard*ldglobalboard * sizeof(int));
    num_alive = generate_initial_board( BS, &globboard[1+ldglobalboard] , ldglobalboard );
    output_board( BS, &globboard[1+ldglobalboard], ldglobalboard, 0 );
    fprintf(stderr, "Starting number of living cells = %d\n", num_alive);
  }

  // Matrix block type used by each processes
  MPI_Datatype block2, block;
  MPI_Type_vector(ldboard-2, ldboard-2, ldglobalboard, MPI_INT, &block2);
  MPI_Type_create_resized(block2, 0, sizeof(int), &block);
  MPI_Type_commit(&block);

  // Matrix sub block type used by each processes
  MPI_Datatype sub_block2, sub_block;
  MPI_Type_vector(ldboard-2, ldboard-2, ldboard, MPI_INT, &sub_block2);
  MPI_Type_create_resized(sub_block2, 0, sizeof(int), &sub_block);
  MPI_Type_commit(&sub_block);

  int *process_count = (int*)malloc(size*sizeof(int));  
  // number of cells per processes
  int *cell_per_processes = (int*)malloc(size*sizeof(int));

  // Prototyping moves for each processes (preparing matrix's scatter)
  for (i = 0; i < process_per_row; ++i){
    for (j = 0; j < process_per_column; ++j){
      process_count[i+j*process_per_column]= 1;
      cell_per_processes[i+j*process_per_column]= i*ldglobalboard*(ldboard-2)+j*(ldboard-2);
    }
  }

  /* Explodes matrix into sub_blocks elements */
  MPI_Scatterv(&globboard[1+ldglobalboard], process_count, cell_per_processes, block, &board[ldboard+1], 1, sub_block,0, grid);

  // Initialize for each processes, a table of the neighbours.
  int neighbours[8];
  neighbour_table(neighbours, grid, rank);

  /* Time to begin */
  t1 = mytimer();
  int blocksize = ldboard-2;
  MPI_Datatype row_blocks;
  MPI_Type_vector(blocksize, 1, ldboard, MPI_INT, &row_blocks);
  MPI_Type_commit(&row_blocks);

  // status for waiting time...
  MPI_Status mpi_status;

  // Create as much MPI request as number of neighbours possible (in the worst case 8) 
  MPI_Request cart_request[8];
  for (loop = 1; loop <= loop_iterations; ++loop) {
    /* Start communications to send and recv informations from neighbours */
    inter_proc_communications(cart_request, neighbours, grid, blocksize, board, ldboard, row_blocks);

    /* Compute inside process cells */
    for (j = 2; j <= blocksize-1; ++j) {
      for (i = 2; i <= blocksize-1; ++i) {
	ngb( i, j ) =
	  cell( i-1, j-1 ) + cell( i, j-1 ) + cell( i+1, j-1 ) +
	  cell( i-1, j   ) +                  cell( i+1, j   ) +
	  cell( i-1, j+1 ) + cell( i, j+1 ) + cell( i+1, j+1 );
      }
    }

    /* Computes cells on the border */

    // Cell neighbour's composition
    
    // 4 2 5       4           4 2 5       4 2 5       4 2 5 //
    // 0 X 1  -->  0      -->  0      -->  0   1  -->  0   1 //
    // 6 3 7       6           6           6   7       6 3 7 //
    
    /* Column on the left needs data from the left process --> 4, 0, 6*/ 
    MPI_Wait(&cart_request[0], &mpi_status);
    MPI_Wait(&cart_request[4], &mpi_status);
    MPI_Wait(&cart_request[6], &mpi_status);
    process_frontier(1, blocksize, board, COLUMN, ldboard, nbngb, ldnbngb);

    /* Line above needs data from the above process --> 2, 5 */
    MPI_Wait(&cart_request[2], &mpi_status);
    MPI_Wait(&cart_request[5], &mpi_status);
    process_frontier(1, blocksize, board, ROW, ldboard, nbngb, ldnbngb);

    /* Column on the right needs data from the right process --> 1, 7 */
    MPI_Wait(&cart_request[1], &mpi_status);
    MPI_Wait(&cart_request[7], &mpi_status);
    process_frontier(blocksize, blocksize, board, COLUMN, ldboard, nbngb, ldnbngb);

    /* Line under needs data from under process --> 3 */
    MPI_Wait(&cart_request[3], &mpi_status);
    process_frontier(blocksize, blocksize, board, ROW, ldboard, nbngb, ldnbngb);


    /* Update the cell */
    num_alive = 0;
    for (j = 1; j <= blocksize; ++j) {
      for (i = 1; i <= blocksize; ++i) {
	if ( (ngb( i, j ) < 2) ||
	     (ngb( i, j ) > 3) ) {
	  cell(i, j) = 0;
	}
	else {
	  if ((ngb( i, j )) == 3)
	    cell(i, j) = 1;
	}
	if (cell(i, j) == 1) {
	  num_alive+=1;
	}
      }
    }
    printf("%d \n", num_alive);
  }

  /* Reassembles matrix into one from the sub blocks in the block */
  MPI_Gatherv(&board[ldboard+1], 1, sub_block, &globboard2[1+ldglobalboard], process_count, cell_per_processes, block, 0, grid);

  /* Reduction to determine max time execution */
  time = mytimer() - t1;
  MPI_Allreduce(&time, &final_time, 1,MPI_DOUBLE, MPI_MAX, grid);
  
  /* Reduction to determine number of cells still alive in all processes */
  MPI_Allreduce(MPI_IN_PLACE, &num_alive, 1, MPI_INT, MPI_SUM, grid);
  
  /* The END */
  if(!rank){
    // Combien de cellules sont en PLS à la fin de la soirée ?
    printf("Final number of living cells = %d\n", num_alive);
    printf("time=%.2lf ms\n",(double)time * 1.e3);
    char fname [40];
    sprintf(fname, "mpi_debug_%d.dat", size);
    FILE* f=fopen(fname, "w");
    // JUST TELL ME IF IT WORKS !!
    if (f != NULL)
      fprintf(f,"%.2lf", time*1.e3);
    fclose(f);
    output_board( BS, &globboard2[1+ldglobalboard], ldglobalboard, loop_iterations);
  }
  // FREE ALL
  free(process_count);
  free(cell_per_processes);
  free(board);
  free(nbngb);
  MPI_Finalize();
  // The final end
  return EXIT_SUCCESS;
}
