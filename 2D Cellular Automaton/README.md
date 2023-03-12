### Files introduction

The code file contains a header file and some source files.

The header file is "automaton.h", which contains the declaration of all the functions.



The source files are "automaton.c", "basic_function.c", "mpi_function.c", cellio.c", "unirand.c". More details will be delivered below:

- automaton.c: This is the main file which should be executed.
- basic_function.c: This is the file which contains most of the basic functions without MPI libraries.
  - set_parameters: set the arguments by terminal inputs.
  - divide_size: automatically create a virtual topology and set the process coordination.
  - init_original_cell: initialize the big cell\[L][L] with random live cells.
  - update_live_cells: update the cells with halo and calculate the number of live cells.
  
- mpi_function.c: This is the file which contains mpi communication functions.
  - set_topology: automatically set the virtual topology and its related variables.
  - divide_cells: divide the big cell\[L][L] to each process
  - send_halos: send halos to neighbour processes.
  - gather_cells: gather cells to rank 0.

- unirand.c: This is the file which is used to generate random numbers.
- cellio.c: This is the file which is used to write down the live cells image.





### How to build the code on Cirrus?

Using the command "make" can compile the code. The code will be build with Intel compilers with -O3, mpt/2.22, and intel-compilers-19 modules.



### How to run the code?

The easiest way to run the code in the command line is type `mpirun -n process_number ./automaton seed_number`. The process_number you want to use. The seed_number is the seed which generate random live cells in the original cell\[L][L]



For additional arguments, users can type `mpirun -n process_number ./automaton` and plus argument type:

- -s seed_number: this is the same as the easiest way which shows above.
- -d density_number: this is to set the density number.
- -l L_number: this is to set the L of the cell\[L][L]