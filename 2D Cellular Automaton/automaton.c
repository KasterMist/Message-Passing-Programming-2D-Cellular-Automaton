#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#include "automaton.h"

int main(int argc, char *argv[])
{
  struct BasicParameter basic;
  struct Neighbour neighbour;

  basic.seed = 1234;
  basic.rho = 0.49;
  basic.L = 768;

  /*
   * Set MPI variables
   */ 
  MPI_Comm comm = MPI_COMM_WORLD;
  int size, rank;
  MPI_Init(&argc, &argv);

  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  /*
   * Deal with argument input in the terminal
   */ 
  set_parameters(argc, argv, &basic, rank); 

  if(rank == 0){
    fprintf(stdout, "seed is %d, density is %f\n", basic.seed, basic.rho);
  }
  
  /*
   * Create the cart topology and define the neighbours
   */
  int dims[2] = {0, 0}, periods[2] = {0, 1}, coords[2], size_coords[4] = {0, 0, 0, 0};
  int ndims = 2, reorder = 0;
  set_topology(size, rank, &comm, dims, size_coords, &neighbour, &basic);

  int allcell[basic.L][basic.L];
  int cell[basic.LX+2][basic.LY+2];

  int ncell, localncell, step, maxstep, printfreq;
  maxstep = 10 * basic.L;
  printfreq = 500;
  int origin_live_cell = 0;
  
  /*
   *  Update for a fixed number of steps, periodically report progress, and separate original cells to each process
   */
  if (rank == 0)
  {
      printf("automaton: running on %d process(es)\n", size);
      printf("automaton: L = %d, rho = %f, seed = %d, maxstep = %d\n",
             basic.L, basic.rho, basic.seed, maxstep);
      
      rinit(basic.seed);

      init_original_cell(basic.rho, &origin_live_cell, basic.L, allcell); 
      printf("automaton: rho = %f, live cells = %d, actual density = %f\n",
          basic.rho, origin_live_cell, ((double) origin_live_cell)/((double) basic.L*basic.L) );
  }
  MPI_Bcast(&origin_live_cell, 1, MPI_INT, 0, comm);
  divide_cells(rank, size, dims, basic, size_coords, allcell, cell, comm);
  

  /*
   * Start the timer here.
   * Calculate and update live cells for each step
   */
  double time;
  time = MPI_Wtime();
  int terminal_flag = 0;

  for (step = 1; step <= maxstep; step++){
    send_halos(basic.LX, basic.LY, cell, neighbour, comm); 
    
    /*
     * Update cells in each process
     */
    localncell = update_live_cells(basic.LX, basic.LY, cell);

    /*
     *  Check whether the number of liver cells reaches the terminal cnodition
     */
    MPI_Allreduce(&localncell, &ncell, 1, MPI_INT, MPI_SUM, comm); 
    if((ncell <= 2*origin_live_cell/3) || (ncell >= 3*origin_live_cell/2)){
      if(rank == 0){
        printf("The code reaches the termination condition. The live cells on step %d is %d\n", step, ncell);
      } 
      break; 
    }
    
    /*
     * Print steps with live cells
     */
    if (step % printfreq == 0){
      if (rank == 0){
        printf("automaton: number of live cells on step %d is %d\n", step, ncell);
      }
    }

  }

  /*
   * Stop the timer
   */
  time = MPI_Wtime() - time;
  if(rank == 0){
    printf("The time is %.4f seconds.\n", time);
  }

  /*
   *  Gather cells back to allcell
   */
  gather_cells(rank, size, dims, basic, size_coords, allcell, cell, comm);

  if (rank == 0){
      cellwrite("cell.pbm", basic.L, allcell);
      printf("wrtie done!\n");
  }

  MPI_Finalize();

  return 0;
}
