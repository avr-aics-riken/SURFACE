#ifdef LSGL_ENABLE_MPI
#include <mpi.h>
#endif

#include <stdio.h>

#include "gtest/gtest.h"

GTEST_API_ int main(int argc, char **argv) {
#ifdef LSGL_ENABLE_MPI
  MPI_Init(&argc, &argv);
#endif
  printf("Running main() from gtest_main.cc\n");
  testing::InitGoogleTest(&argc, argv);
  int ret = RUN_ALL_TESTS();
#ifdef LSGL_ENABLE_MPI
  MPI_Finalize();
#endif
  return ret;
}
