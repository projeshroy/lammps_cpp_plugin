#pragma once

#include <string>
#include <fstream>
#include <iostream>
#include <ios>
#include <iomanip>
#include <streambuf>
#include <map>
#include <set>
#include <array>
//#include <complex.h>
#include <random>
#include <type_traits>
#include <exception>
#include <algorithm>
#include <chrono>
#include <thread>
#include <atomic>
#include <mutex>
#include <tuple>
#include <omp.h>
#include <time.h>
#include <cmath>
#include <cassert>
#include <math.h>

#define EIGEN_DEFAULT_TO_COLUMN_MAJOR
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues> 
#include "EigenDummyTensor.h"
#include "EigenMultiDimArray.h"

#define LAMMPS_LIB_MPI
#include "/home/projesh/lammps/src/lammps.h"
#include "/home/projesh/lammps/src/library.h"
#include "/home/projesh/lammps/src/input.h"
#include "/home/projesh/lammps/src/atom.h"
#include "/home/projesh/openmpi-4.1.5/build_gpu/include/mpi.h"
using namespace LAMMPS_NS;


