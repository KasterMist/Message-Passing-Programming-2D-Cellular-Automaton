#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "automaton.h"

void set_topology(int size, int rank, MPI_Comm *comm, int dims[2], int size_coords[4],
 struct Neighbour *neighbour, struct BasicParameter *basic){
  
  int ndims = 2, reorder = 0;
  int periods[2] = {0, 1};
  int coords[2];

  MPI_Dims_create(size, ndims, dims);
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, comm);
  
  
  MPI_Cart_coords(*comm, rank, 2, coords);
  MPI_Cart_shift(*comm, 0, 1, &(*neighbour).left, &(*neighbour).right);    
  MPI_Cart_shift(*comm, 1, 1, &(*neighbour).down, &(*neighbour).up);

  
  divide_size((*basic).L, dims, coords, size_coords);

  if(coords[0] == dims[0] - 1){
    (*basic).LX = size_coords[2];
  }else{
    (*basic).LX = size_coords[0]; 
  }

  if(coords[1] == dims[1] - 1){
    (*basic).LY = size_coords[3]; 
  }else{
    (*basic).LY = size_coords[1];
  }

}


void divide_cells(int rank, int size, int dims[2], struct BasicParameter basic,
 int size_coords[4], int allcell[basic.L][basic.L], int cell[basic.LX+2][basic.LY+2], MPI_Comm comm){
  
  int L = basic.L;
  int LX = basic.LX;
  int LY = basic.LY;

  MPI_Datatype normal_vector_type;
  MPI_Type_vector(size_coords[0], size_coords[1], L, MPI_INT, &normal_vector_type);
  MPI_Type_commit(&normal_vector_type);
  
  MPI_Datatype long_x_vector_type;
  MPI_Type_vector(size_coords[2], size_coords[1], L, MPI_INT, &long_x_vector_type);
  MPI_Type_commit(&long_x_vector_type);

  MPI_Datatype long_y_vector_type;
  MPI_Type_vector(size_coords[0], size_coords[3], L, MPI_INT, &long_y_vector_type);
  MPI_Type_commit(&long_y_vector_type); 

  MPI_Datatype special_vector_type;
  MPI_Type_vector(size_coords[2], size_coords[3], L, MPI_INT, &special_vector_type);
  MPI_Type_commit(&special_vector_type);
   
  MPI_Datatype small_cell_vector_type;
  MPI_Type_vector(LX, LY, LY+2, MPI_INT, &small_cell_vector_type);
  MPI_Type_commit(&small_cell_vector_type); 

  MPI_Status recv_status; 
  int recv_coords[2];
  if(rank == 0){
    for(int i = 0; i < LX; i++){
      for(int j = 0; j < LY; j++){
        cell[i+1][j+1] = allcell[i][j];
      }
    }
    for(int sr = 1; sr < size; sr++){
      MPI_Cart_coords(comm, sr, 2, recv_coords);

      if(sr == size - 1){
        MPI_Ssend(&allcell[LX*recv_coords[0]][LY*recv_coords[1]], 1, special_vector_type, sr, 0, comm); 
      }
      else if(recv_coords[0] == dims[0] - 1){
        MPI_Ssend(&allcell[LX*recv_coords[0]][LY*recv_coords[1]], 1, long_x_vector_type, sr, 0, comm);
      }
      else if(recv_coords[1] == dims[1] - 1){
        MPI_Ssend(&allcell[LX*recv_coords[0]][LY*recv_coords[1]], 1, long_y_vector_type, sr, 0, comm); 
      }
      else{
        MPI_Ssend(&allcell[LX*recv_coords[0]][LY*recv_coords[1]], 1, normal_vector_type, sr, 0, comm);
      }
      
    } 
  }
  else{
    MPI_Recv(&cell[1][1], 1, small_cell_vector_type, 0, 0, comm, &recv_status);
  }
}


void send_halos(int LX, int LY, int cell[LX+2][LY+2], struct Neighbour neighbour, MPI_Comm comm){

  MPI_Datatype boundary_vector_type;
  MPI_Type_vector(LX, 1, LY+2, MPI_INT, &boundary_vector_type);
  MPI_Type_commit(&boundary_vector_type);
  
  MPI_Request send_request_list[4];
  MPI_Request recv_request_list[4];
  MPI_Status send_status_list[4];
  MPI_Status recv_status_list[4];
  MPI_Status send_status;
  MPI_Status recv_status;

  int send_left_tag = 1, send_right_tag = 2, send_down_tag = 3, send_up_tag = 4;

  MPI_Isend(&cell[1][1], LY, MPI_INT, neighbour.left, send_left_tag, comm, &send_request_list[0]);      
  MPI_Isend(&cell[LX][1], LY, MPI_INT, neighbour.right, send_right_tag, comm, &send_request_list[1]);
  MPI_Isend(&cell[1][1], 1, boundary_vector_type, neighbour.down, send_down_tag, comm, &send_request_list[2]);
  MPI_Isend(&cell[1][LY], 1, boundary_vector_type, neighbour.up, send_up_tag, comm, &send_request_list[3]); 

  MPI_Irecv(&cell[0][1], LY, MPI_INT, neighbour.left, send_right_tag, comm, &recv_request_list[0]);
  MPI_Irecv(&cell[LX+1][1], LY, MPI_INT, neighbour.right, send_left_tag, comm, &recv_request_list[1]);
  MPI_Irecv(&cell[1][0], 1, boundary_vector_type, neighbour.down, send_up_tag, comm, &recv_request_list[2]); 
  MPI_Irecv(&cell[1][LY+1], 1, boundary_vector_type, neighbour.up, send_down_tag, comm, &recv_request_list[3]);

  MPI_Waitall(4, send_request_list, send_status_list);
  MPI_Waitall(4, recv_request_list, recv_status_list);

}


void gather_cells(int rank, int size, int dims[2], struct BasicParameter basic, 
 int size_coords[4], int allcell[basic.L][basic.L], int cell[basic.LX+2][basic.LY+2], MPI_Comm comm){

  int L = basic.L;
  int LX = basic.LX;
  int LY = basic.LY;
  MPI_Status status;

  MPI_Datatype normal_vector_type;
  MPI_Type_vector(size_coords[0], size_coords[1], L, MPI_INT, &normal_vector_type);
  MPI_Type_commit(&normal_vector_type);
  
  MPI_Datatype long_x_vector_type;
  MPI_Type_vector(size_coords[2], size_coords[1], L, MPI_INT, &long_x_vector_type);
  MPI_Type_commit(&long_x_vector_type);

  MPI_Datatype long_y_vector_type;
  MPI_Type_vector(size_coords[0], size_coords[3], L, MPI_INT, &long_y_vector_type);
  MPI_Type_commit(&long_y_vector_type); 

  MPI_Datatype special_vector_type;
  MPI_Type_vector(size_coords[2], size_coords[3], L, MPI_INT, &special_vector_type);
  MPI_Type_commit(&special_vector_type);
   
  MPI_Datatype small_cell_vector_type;
  MPI_Type_vector(LX, LY, LY+2, MPI_INT, &small_cell_vector_type);
  MPI_Type_commit(&small_cell_vector_type); 

  // sr means subject rank
  if(rank != 0){
    MPI_Ssend(&cell[1][1], 1, small_cell_vector_type, 0, 0, comm);
  }else{
    int recv_coords[2];
    for(int tr = 1; tr < size; tr++){
      MPI_Cart_coords(comm, tr, 2, recv_coords);

      if(tr == size - 1){
        MPI_Recv(&allcell[LX*recv_coords[0]][LY*recv_coords[1]], 1, special_vector_type, tr, 0, comm, &status); 
      }
      else if(recv_coords[0] == dims[0] - 1){
        MPI_Recv(&allcell[LX*recv_coords[0]][LY*recv_coords[1]], 1, long_x_vector_type, tr, 0, comm, &status); 
      }
      else if(recv_coords[1] == dims[1] - 1){
        MPI_Recv(&allcell[LX*recv_coords[0]][LY*recv_coords[1]], 1, long_y_vector_type, tr, 0, comm, &status); 
      }
      else{
        MPI_Recv(&allcell[LX*recv_coords[0]][LY*recv_coords[1]], 1, normal_vector_type, tr, 0, comm, &status); 
      }
      
    }

    for(int i = 0; i < LX; i++){
      for(int j = 0; j < LY; j++){
        allcell[i][j] = cell[i+1][j+1];
      }
    }
    
  }
}