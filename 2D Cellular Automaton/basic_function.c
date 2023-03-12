#include <stdio.h>
#include <stdlib.h>

#include "automaton.h"

void set_parameters(int argc, char *argv[], struct BasicParameter *basic, int rank){
  for(int i = 1; i < argc; i++){
    switch (argv[i][1])
    {
    case 's':
      if(i+1 == argc) {
        if(rank == 0) printf("Input format error!\n"); 
        MPI_Finalize();
        exit(1);
      }
      (*basic).seed = atoi(argv[i+1]); i++; break;
    case 'd':
      if(i+1 == argc){
        if(rank == 0) printf("Input format error!\n");
        MPI_Finalize();; 
        exit(1);
      }
      (*basic).rho = atof(argv[i+1]); i++; break;
    case 'l':
      if(i+1 == argc){ 
        if(rank == 0) printf("Input format error!\n");
        MPI_Finalize();
        exit(1);
      }
      (*basic).L = atoi(argv[i+1]); i++; break;
    default:
      if(argc > 2){
        if(rank == 0) printf("Input format error!\n");
        MPI_Finalize();
        exit(1);
      }
      (*basic).seed = atoi(argv[i]);
      break;
    }
    
  }
}


void divide_size(int L, int dims[2], int coords[2], int size_coords[4]){
  int x_quotient = L / dims[0];
  int x_remainder = L % dims[0];
  int y_quotient = L / dims[1];
  int y_remainder = L % dims[1];
  
  size_coords[0] = x_quotient;
  size_coords[1] = y_quotient;
  size_coords[2] = x_quotient + x_remainder;
  size_coords[3] = y_quotient + y_remainder;
  
}


void init_original_cell(double rho, int *ncell, int l, int output_cell[l][l]){
  
  double r;

  int n = 0;

  for (int i=0; i < l; i++){
    for (int j=0; j < l; j++){
      r = uni();

      if(r < rho){
        output_cell[i][j] = 1;
        n++;
      }
      else{
        output_cell[i][j] = 0;
      }
    }
  }

 *ncell = n; 
}


int update_live_cells(int LX, int LY, int cell[LX+2][LY+2]){
  int neigh[LX+2][LY+2];
  for (int i = 1; i <= LX; i++){
      for (int j = 1; j <= LY; j++){
         /*
          * Set neigh[i][j] to be the sum of cell[i][j] plus its
          * four nearest neighbours
          */

        neigh[i][j] =   cell[i][j] 
                      + cell[i][j-1]
                      + cell[i][j+1]
                      + cell[i-1][j]
                      + cell[i+1][j];
      }
    }

    int localncell = 0;

    for (int i = 1; i <= LX; i++){
      for (int j = 1; j <= LY; j++){
        /*
          * Update based on number of neighbours
          */

        if (neigh[i][j] == 2 || neigh[i][j] == 4 || neigh[i][j] == 5){
          cell[i][j] = 1;
          localncell++;
        }
        else{
          cell[i][j] = 0;
        }
      }
    }
    return localncell;
}