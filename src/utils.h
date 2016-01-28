#ifndef __UTIL_H__
#define __UTIL_H__

#include <mpi.h>

typedef struct cart_grid{
    int rank_i;      // rank row
    int rank_j;      // rank column
    int above_p;     // (i+1,j) 
    int under_p;     // (i-1,j)
    int left_p;      // (i,j-1)
    int right_p;     // (i,j+1)
} cart_grid;

enum dispatch{GATHER, SCATTER};

int init_cart_grid(int rank, MPI_Comm *communicator, cart_grid *grid);

int dispatch_matrix(MPI_Comm *communicator, int *sbuffer, int* rbuffer, int row_procs_nb, int blocksize, enum dispatch kind, int block);

int make_communicator(int tot_procs, int *row_procs, MPI_Comm *communicator, int *rank);

#endif // __UTIL_H__
