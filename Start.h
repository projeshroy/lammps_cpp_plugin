#pragma once
#include "declarations.h"
#include "includeallglobals.h"

void Start(int argc, char **argv){

	lmp = new LAMMPS(argc, (char **)argv, MPI_COMM_WORLD);

	if(mpi_id == 0)
	std::cout << " WARNING! There are some issues with GPU packages. Use OMP instead" << std::endl;

	if (lmp == NULL) {
	printf("LAMMPS initialization failed");
	throw std::exception();
  	}

}

