#include "gtest/gtest.h"

TEST(DummyTest, Dummy) {
  double a = 0.1;
  double b = 0.1;

  EXPECT_DOUBLE_EQ(a, b);
}
