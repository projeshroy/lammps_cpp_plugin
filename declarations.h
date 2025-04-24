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
#include <lammps.h>
#include <library.h>
#include <input.h>
#include <atom.h>
#include <mpi.h>

using namespace LAMMPS_NS;

void omp_set_num_threads(int num_threads);

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec_d;
typedef Eigen::Matrix<int, Eigen::Dynamic, 1> Vec_i;
typedef Eigen::Matrix<unsigned int, Eigen::Dynamic, 1> Vec_ui;
typedef Eigen::Matrix<std::string, Eigen::Dynamic, 1> Vec_s;
typedef Eigen::Matrix<bool, Eigen::Dynamic, 1> Vec_b;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> Vec_Compl_d;

typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> Mat_i;
typedef Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic> Mat_ui;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Mat_d;
typedef Eigen::Matrix<std::string, Eigen::Dynamic, Eigen::Dynamic> Mat_s;
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> Mat_b;
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> Mat_Compl_d;
