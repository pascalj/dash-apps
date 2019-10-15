#pragma once

#include <mpi.h>
#include <iomanip>

double t(double start, const char* msg) {
  double end = MPI_Wtime();
  std::cout << "[" << std::setprecision(4) << end - start << "] " << msg << std::endl;
  return end;
}
