#include "includeallglobals.h"

void Start(int argc, char **argv){
	MPI_Init(&argc, &argv);
	LAMMPS_HANDLE = lammps_open(0, NULL, MPI_COMM_WORLD, NULL);
	
	if (LAMMPS_HANDLE == NULL) {
	printf("LAMMPS initialization failed");
	throw std::exception();
  	}
}

