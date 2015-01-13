#include "gtest/gtest.h"

#include "matrix.h"

using namespace lsgl::render;

TEST(MatrixTest, Identity) {
  double a = 0.1;
  double b = 0.1;

  double m[4][4];

  Matrixd::Identity(m);

  EXPECT_DOUBLE_EQ(1.0, m[0][0]);
  EXPECT_DOUBLE_EQ(1.0, m[1][1]);
  EXPECT_DOUBLE_EQ(1.0, m[2][2]);
  EXPECT_DOUBLE_EQ(1.0, m[3][3]);

  EXPECT_DOUBLE_EQ(0.0, m[0][1]);
}
