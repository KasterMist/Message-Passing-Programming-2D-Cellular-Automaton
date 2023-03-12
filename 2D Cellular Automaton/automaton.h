// #ifndef MPI_INCLUDED
// #define MPI_INCLUDED
#include "mpi.h"
/*
 *  Main header file for percolation exercise.
 */

struct Neighbour{
  int left;
  int right;
  int down;
  int up;
};

struct BasicParameter{
  int seed;
  double rho;
  int L;
  int LX;
  int LY;
};

void cellwrite(char *cellfile, int l, int cell[l][l]);
void cellwritedynamic(char *cellfile, int **cell, int l);

/*
 *  Random numbers
 */
void rinit(int ijkl);
float uni(void);

/*
 *  Basic functions
 */
void set_parameters(int argc, char *argv[], struct BasicParameter *basic, int rank);
void divide_size(int L, int dims[2], int coords[2], int size_coords[4]);
void init_original_cell(double rho, int *ncell, int l, int output_cell[l][l]);
int update_live_cells(int LX, int LY, int cell[LX+2][LY+2]);


/*
 *  MPI functions
 */
void set_topology(int size, int rank, MPI_Comm *comm, int dims[2], int size_coords[4],
 struct Neighbour *neighbour, struct BasicParameter *basic);

void divide_cells(int rank, int size, int dims[2], struct BasicParameter basic,
  int size_coords[4], int allcell[basic.L][basic.L], int cell[basic.LX+2][basic.LY+2], MPI_Comm comm);

void send_halos(int LX, int LY, int cell[LX+2][LY+2], struct Neighbour neighbour, MPI_Comm comm);

void gather_cells(int rank, int size, int dims[2], struct BasicParameter basic, 
  int size_coords[4], int allcell[basic.L][basic.L], int cell[basic.LX+2][basic.LY+2], MPI_Comm comm);


// #endif