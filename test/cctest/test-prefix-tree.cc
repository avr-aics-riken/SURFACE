#include "gtest/gtest.h"

#include "../render/render_common.h"
#include "../render/tinymt64.h"
#include "../render/prefix_tree_util.h"

using namespace lsgl::render;

namespace {
tinymt64_t rng;

void InitRandom(int seed) {
  rng.mat1 = 0xfa051f40; 
  rng.mat2 = 0xffd0fff4;
  rng.tmat = 0x58d02ffeffbfffbcULL;
  tinymt64_init(&rng, seed);
}

float GenRandomReal()
{
  return (float)(tinymt64_generate_double(&rng));

}

} // namespace

TEST(LSGLRenderTest, TestMortonCode30) {

  assert(sizeof(real) == sizeof(float));

  real3 one_point_5(0.5, 0.5, 0.5);
  real3 bmin(0.0, 0.0, 0.0);
  real3 bmax(1.0, 1.0, 1.0);
  real  inv = 1024.0;

  uint32_t code = MortionCode30(one_point_5, bmin, inv, inv, inv);

  uint32_t expect = (1<<29) | (1<<28) | (1<<27);
  EXPECT_EQ(expect, code);

}

#if defined(__SSE2__) || (_M_IX86_FP >= 2)
TEST(LSGLRenderTest, TestComputeMortionCodesSIMD) {
  assert(sizeof(real) == sizeof(float));

  int ns[] = {8, 16, 32};

  int niter = sizeof(ns) / sizeof(int);

  int seed = 123456;
  InitRandom(seed);

  for (int i = 0; i < niter; i++) {

    int npoints = ns[i];

    std::vector<float> points(3*npoints);
    std::vector<int> indices(npoints);
    std::vector<uint32_t> codes(npoints);
    std::vector<uint32_t> codes2(npoints);

    for (int j = 0; j < npoints; j++) {
      points[3*j+0] = GenRandomReal();
      points[3*j+1] = GenRandomReal();
      points[3*j+2] = GenRandomReal();

      indices[j] = j;
    }

    real3 bmin(0.0, 0.0, 0.0);
    real3 bmax(1.0, 1.0, 1.0);

    CalculateMortonCodes30(&codes.at(0), &points.at(0), bmin, bmax, 0, npoints);
    CalculateMortonCodes30SIMD(&codes2.at(0), &points.at(0), bmin, bmax, 0, npoints);

    for (int j = 0; j < npoints; j++) {
      EXPECT_EQ(codes[j], codes2[j]);
    }

  }
}
#endif


TEST(LSGLRenderTest, TestRadixSort) {

  assert(sizeof(real) == sizeof(float));

  int ns[] = {8, 16, 32};

  int niter = sizeof(ns) / sizeof(int);

  int seed = 123456;
  InitRandom(seed);

  for (int i = 0; i < niter; i++) {

    int npoints = ns[i];

    std::vector<float> points(3*npoints);
    std::vector<int> indices(npoints);
    std::vector<uint32_t> codes(npoints);

    for (int j = 0; j < npoints; j++) {
      points[3*j+0] = GenRandomReal();
      points[3*j+1] = GenRandomReal();
      points[3*j+2] = GenRandomReal();

      indices[j] = j;
    }

    real3 bmin(0.0, 0.0, 0.0);
    real3 bmax(1.0, 1.0, 1.0);

    CalculateMortonCodes30(&codes.at(0), &points.at(0), bmin, bmax, 0, npoints);

    std::vector<IndexKey30> keys(npoints);
    for (size_t i = 0; i < npoints; i++) {
      keys[i].index = indices[i];
      keys[i].code  = codes[i];
    }
    RadixSortByMortionCode30MSB(&keys.at(0), &keys.at(0) + npoints);

    for (size_t i = 0; i < npoints-1; i++) {
      EXPECT_LE(keys[i].code, keys[i+1].code); 
      if (keys[i].code == keys[i+1].code) {
        EXPECT_LE(keys[i].index, keys[i+1].index); 
      }
    }
  }
}

TEST(LSGLRenderTest, TestRadixSort30) {

  assert(sizeof(real) == sizeof(float));

  int ns[] = {8, 16, 32, 1023, 11142, 3456663};

  int niter = sizeof(ns) / sizeof(int);

  int seed = 123456;
  InitRandom(seed);

  for (int i = 0; i < niter; i++) {

    int npoints = ns[i];

    std::vector<float> points(3*npoints);
    std::vector<int> indices(npoints);
    std::vector<uint32_t> codes(npoints);

    for (int j = 0; j < npoints; j++) {
      points[3*j+0] = GenRandomReal();
      points[3*j+1] = GenRandomReal();
      points[3*j+2] = GenRandomReal();

      indices[j] = j;
    }

    real3 bmin(0.0, 0.0, 0.0);
    real3 bmax(1.0, 1.0, 1.0);

    CalculateMortonCodes30(&codes.at(0), &points.at(0), bmin, bmax, 0, npoints);

    std::vector<IndexKey30> keys(npoints);
    for (size_t i = 0; i < npoints; i++) {
      keys[i].index = indices[i];
      keys[i].code  = codes[i];
    }

    //std::vector<IndexKey30> temp = keys;
    //MergeSort30(&keys.at(0), &temp.at(0), 0, npoints-1);
    RadixSort30(&keys.at(0), &keys.at(0) + npoints);

    for (size_t i = 0; i < npoints-1; i++) {
      //std::string k = BitString32(keys[i].code);
      //printf("[%010d] = %010d(%s)\n", i, keys[i].code, k.c_str());
      EXPECT_LE(keys[i].code, keys[i+1].code); 
      //if (keys[i].code == keys[i+1].code) {
      //  EXPECT_LE(keys[i].index, keys[i+1].index); 
      //}
    }
  }
}

TEST(LSGLRenderTest, TestParallelRadixSort30) {

  assert(sizeof(real) == sizeof(float));

  int ns[] = {8, 16, 32, 1023, 11142, 3456663};

  int niter = sizeof(ns) / sizeof(int);

  int seed = 123456;
  InitRandom(seed);

  for (int i = 0; i < niter; i++) {

    int npoints = ns[i];

    std::vector<float> points(3*npoints);
    std::vector<int> indices(npoints);
    std::vector<uint32_t> codes(npoints);

    for (int j = 0; j < npoints; j++) {
      points[3*j+0] = GenRandomReal();
      points[3*j+1] = GenRandomReal();
      points[3*j+2] = GenRandomReal();

      indices[j] = j;
    }

    real3 bmin(0.0, 0.0, 0.0);
    real3 bmax(1.0, 1.0, 1.0);

    CalculateMortonCodes30(&codes.at(0), &points.at(0), bmin, bmax, 0, npoints);

    std::vector<IndexKey30> keys(npoints);
    for (size_t i = 0; i < npoints; i++) {
      keys[i].index = indices[i];
      keys[i].code  = codes[i];
    }

    //std::vector<IndexKey30> temp = keys;
    //MergeSort30(&keys.at(0), &temp.at(0), 0, npoints-1);
    RadixSort30(&keys.at(0), &keys.at(0) + npoints);

    for (size_t i = 0; i < npoints-1; i++) {
      //std::string k = BitString32(keys[i].code);
      //printf("[%010d] = %010d(%s)\n", i, keys[i].code, k.c_str());
      EXPECT_LE(keys[i].code, keys[i+1].code); 
      //if (keys[i].code == keys[i+1].code) {
      //  EXPECT_LE(keys[i].index, keys[i+1].index); 
      //}
    }
  }
}

TEST(LSGLRenderTest, TestMergeSort) {

  assert(sizeof(real) == sizeof(float));

  int ns[] = {8, 16, 32, 1023, 11142, 3456663};

  int niter = sizeof(ns) / sizeof(int);

  int seed = 123456;
  InitRandom(seed);

  for (int i = 0; i < niter; i++) {

    int npoints = ns[i];

    std::vector<float> points(3*npoints);
    std::vector<int> indices(npoints);
    std::vector<uint32_t> codes(npoints);

    for (int j = 0; j < npoints; j++) {
      points[3*j+0] = GenRandomReal();
      points[3*j+1] = GenRandomReal();
      points[3*j+2] = GenRandomReal();

      indices[j] = j;
    }

    real3 bmin(0.0, 0.0, 0.0);
    real3 bmax(1.0, 1.0, 1.0);

    CalculateMortonCodes30(&codes.at(0), &points.at(0), bmin, bmax, 0, npoints);

    std::vector<IndexKey30> keys(npoints);
    for (size_t i = 0; i < npoints; i++) {
      keys[i].index = indices[i];
      keys[i].code  = codes[i];
    }

    std::vector<IndexKey30> temp = keys;

    MergeSort30(&keys.at(0), &temp.at(0), 0, npoints-1);

    for (size_t i = 0; i < npoints-1; i++) {
      //std::string k = BitString32(keys[i].code);
      //printf("[%010d] = %010d(%s)\n", i, keys[i].code, k.c_str());
      EXPECT_LE(keys[i].code, keys[i+1].code); 
      //if (keys[i].code == keys[i+1].code) {
      //  EXPECT_LE(keys[i].index, keys[i+1].index); 
      //}
    }
  }
}

TEST(LSGLRenderTest, TestSorter) {

  assert(sizeof(real) == sizeof(float));

  int ns[] = {8, 16, 32, 1023, 11142, 3456663};

  int niter = sizeof(ns) / sizeof(int);

  int seed = 123456;
  InitRandom(seed);

  for (int i = 0; i < niter; i++) {

    int npoints = ns[i];

    std::vector<float> points(3*npoints);
    std::vector<int> indices(npoints);
    std::vector<uint32_t> codes(npoints);

    for (int j = 0; j < npoints; j++) {
      points[3*j+0] = GenRandomReal();
      points[3*j+1] = GenRandomReal();
      points[3*j+2] = GenRandomReal();

      indices[j] = j;
    }

    real3 bmin(0.0, 0.0, 0.0);
    real3 bmax(1.0, 1.0, 1.0);

    CalculateMortonCodes30(&codes.at(0), &points.at(0), bmin, bmax, 0, npoints);

    std::vector<IndexKey30> keys0(npoints);
    for (size_t i = 0; i < npoints; i++) {
      keys0[i].index = indices[i];
      keys0[i].code  = codes[i];
    }

    std::vector<IndexKey30> temp = keys0;

    MergeSort30(&keys0.at(0), &temp.at(0), 0, npoints-1);

    std::vector<IndexKey30> keys1(npoints);
    for (size_t i = 0; i < npoints; i++) {
      keys1[i].index = indices[i];
      keys1[i].code  = codes[i];
    }

    RadixSortByMortionCode30MSB(&keys1.at(0), &keys1.at(0) + npoints);

    for (size_t i = 0; i < npoints; i++) {
      EXPECT_EQ(keys0[i].code, keys1[i].code); 
      //EXPECT_EQ(keys0[i].index, keys1[i].index); 
    }
  }
}

namespace {

void VisitAndCountDepth32(
  int& maxDepth,
  std::vector<NodeInfo32>& nodes,
  unsigned int root,
  int depth)
{
  //printf("root = %d\n", root);
  ASSERT_GT(nodes.size(), root);

  uint32_t leftIndex = nodes[root].childIndex;
  uint32_t rightIndex = nodes[root].childIndex+1;
  bool     isLeftLeaf = (nodes[root].leftType == NODE_TYPE_LEAF) ? true : false;
  bool     isRightLeaf = (nodes[root].rightType == NODE_TYPE_LEAF) ? true : false;

  if (maxDepth <= depth) {
    maxDepth = depth;
  }

  if (isLeftLeaf) {
    // do nothing
  } else {
    ASSERT_GT(nodes.size(), leftIndex);
    VisitAndCountDepth32(maxDepth, nodes, leftIndex, depth+1);
  }

  if (isRightLeaf) {
    // do nothing
  } else {
    ASSERT_GT(nodes.size(), rightIndex);
    VisitAndCountDepth32(maxDepth, nodes, rightIndex, depth+1);
  }
}

void VisitAndCountDepth64(
  int& maxDepth,
  std::vector<NodeInfo64>& nodes,
  unsigned int root,
  int depth)
{
  //printf("root = %d\n", root);
  ASSERT_GT(nodes.size(), root);

  uint32_t leftIndex = nodes[root].index;
  uint32_t rightIndex = nodes[root].index+1;
  bool     isLeftLeaf = (nodes[root].leftType == NODE_TYPE_LEAF) ? true : false;
  bool     isRightLeaf = (nodes[root].rightType == NODE_TYPE_LEAF) ? true : false;

  if (maxDepth <= depth) {
    maxDepth = depth;
  }

  if (isLeftLeaf) {
    // do nothing
  } else {
    ASSERT_GT(nodes.size(), leftIndex);
    VisitAndCountDepth64(maxDepth, nodes, leftIndex, depth+1);
  }

  if (isRightLeaf) {
    // do nothing
  } else {
    ASSERT_GT(nodes.size(), rightIndex);
    VisitAndCountDepth64(maxDepth, nodes, rightIndex, depth+1);
  }
}

void VisitAndCountLeafs32(
  std::vector<int>&      counts,
  std::vector<NodeInfo32>& nodes,
  int root)
{
  //printf("root = %d\n", root);

  ASSERT_LE(0, root);
  ASSERT_GT(nodes.size(), root);

  uint32_t leftIndex = nodes[root].childIndex;
  uint32_t rightIndex = nodes[root].childIndex+1;
  bool     isLeftLeaf = (nodes[root].leftType == NODE_TYPE_LEAF) ? true : false;
  bool     isRightLeaf = (nodes[root].rightType == NODE_TYPE_LEAF) ? true : false;


  //printf("root = %d, leftIdx = %d, leaf = %d / %d\n", root, leftIndex, isLeftLeaf, isRightLeaf);


  if (isLeftLeaf) {
    ASSERT_GT(counts.size(), leftIndex);
    counts[leftIndex]++;
  } else {
    ASSERT_GT(nodes.size(), leftIndex);
    VisitAndCountLeafs32(counts, nodes, leftIndex);
  }

  if (isRightLeaf) {
    ASSERT_GT(counts.size(), rightIndex);
    counts[rightIndex]++;
  } else {
    ASSERT_GT(nodes.size(), rightIndex);
    VisitAndCountLeafs32(counts, nodes, rightIndex);
  }
}

void VisitAndCountLeafs64(
  std::vector<int>&      counts,
  std::vector<NodeInfo64>& nodes,
  int root)
{
  //printf("root = %d\n", root);

  ASSERT_LE(0, root);
  ASSERT_GT(nodes.size(), root);

  uint32_t leftIndex = nodes[root].index;
  uint32_t rightIndex = nodes[root].index+1;
  bool     isLeftLeaf = (nodes[root].leftType == NODE_TYPE_LEAF) ? true : false;
  bool     isRightLeaf = (nodes[root].rightType == NODE_TYPE_LEAF) ? true : false;


  //printf("root = %d, leftIdx = %d, leaf = %d / %d\n", root, leftIndex, isLeftLeaf, isRightLeaf);


  if (isLeftLeaf) {
    ASSERT_GT(counts.size(), leftIndex);
    counts[leftIndex]++;
  } else {
    ASSERT_GT(nodes.size(), leftIndex);
    VisitAndCountLeafs64(counts, nodes, leftIndex);
  }

  if (isRightLeaf) {
    ASSERT_GT(counts.size(), rightIndex);
    counts[rightIndex]++;
  } else {
    ASSERT_GT(nodes.size(), rightIndex);
    VisitAndCountLeafs64(counts, nodes, rightIndex);
  }
}

void CountLeafs32(
  std::vector<int>&      counts,
  std::vector<NodeInfo32>& nodes)
{

  ASSERT_EQ(counts.size()-1, nodes.size());

  for (size_t i = 0; i < nodes.size(); i++) {
    uint32_t index = nodes[i].childIndex;
    bool     isLeftLeaf = (nodes[i].leftType == NODE_TYPE_LEAF) ? true : false;
    bool     isRightLeaf = (nodes[i].rightType == NODE_TYPE_LEAF) ? true : false;
    if (isLeftLeaf) {
      counts[index]++;
    }

    if (isRightLeaf) {
      counts[index+1]++;
    }

  }
}

void CountLeafs64(
  std::vector<int>&      counts,
  std::vector<NodeInfo64>& nodes)
{

  ASSERT_EQ(counts.size()-1, nodes.size());

  for (size_t i = 0; i < nodes.size(); i++) {
    uint32_t index = nodes[i].index;
    bool     isLeftLeaf = (nodes[i].leftType == NODE_TYPE_LEAF) ? true : false;
    bool     isRightLeaf = (nodes[i].rightType == NODE_TYPE_LEAF) ? true : false;
    if (isLeftLeaf) {
      counts[index]++;
    }

    if (isRightLeaf) {
      counts[index+1]++;
    }

  }
}

}

TEST(LSGLRenderTest, TestBinaryPrefixTree30) {

  assert(sizeof(real) == sizeof(float));

  //int ns[] = {8, 16, 32, 1024, 1024*1024, 5555555, 123456789};
  //int ns[] = {8, 15, 31, 1023, 1000, 3211,  4000000, 12345678};
  int ns[] = {8, 15, 31, 1023, 1000, 3211};
  //int ns[] = {8, 16, 31};

  int niter = sizeof(ns) / sizeof(int);

  int seed = 123456;
  InitRandom(seed);

  for (int i = 0; i < niter; i++) {

    int npoints = ns[i];

    std::vector<float> points(3*npoints);
    std::vector<int> indices(npoints);
    std::vector<uint32_t> codes(npoints);

    for (int j = 0; j < npoints; j++) {
      points[3*j+0] = GenRandomReal();
      points[3*j+1] = GenRandomReal();
      points[3*j+2] = GenRandomReal();

      indices[j] = j;
    }

    real3 bmin(0.0, 0.0, 0.0);
    real3 bmax(1.0, 1.0, 1.0);

    CalculateMortonCodes30(&codes.at(0), &points.at(0), bmin, bmax, 0, npoints);

    std::vector<IndexKey30> keys(npoints);
    for (size_t i = 0; i < npoints; i++) {
      keys[i].index = indices[i];
      keys[i].code  = codes[i];
    }
    RadixSortByMortionCode30MSB(&keys.at(0), &keys.at(0) + npoints);

    // Tree size = npoints - 1
    std::vector<NodeInfo32> nodes(npoints-1);
    for (size_t i = 0; i < npoints - 1; i++) {
      nodes[i] = ConstructBinaryRadixTree30(&keys.at(0), i, npoints);
    }

    // Tree must visit each leaf node only once.
    {
      std::vector<int> counter(npoints);
      for (size_t i = 0; i < counter.size(); i++) {
        counter[i] = 0;
      }
    
      CountLeafs32(counter, nodes);

      //for (size_t i = 0; i < counter.size(); i++) {
      //    printf("count[%d] = %d\n", i, counter[i]);
      //}

      for (size_t i = 0; i < counter.size(); i++) {
        if (counter[i] != 1) {
          printf("i = %d, node = %d, leaf = %d, %d\n", i, nodes[i].childIndex, nodes[i].leftType, nodes[i].rightType);
        }
        EXPECT_EQ(1, counter[i]);
      }
    }

    // Leaf tree must be accessible by traversing 
    {
      std::vector<int> counter(npoints);
      for (size_t i = 0; i < counter.size(); i++) {
        counter[i] = 0;
      }
    
      VisitAndCountLeafs32(counter, nodes, 0);

      for (size_t i = 0; i < counter.size(); i++) {
        if (counter[i] != 1) {
          printf("i = %d, node = %d, leaf = %d, %d\n", i, nodes[i].childIndex, nodes[i].leftType, nodes[i].rightType);
        }
        EXPECT_EQ(1, counter[i]);
      }
    }

    // Check depth
    int maxDepth = 0;
    VisitAndCountDepth32(maxDepth, nodes, 0, 0);
    EXPECT_GT(30, maxDepth);
    //printf("maxDepth = %d\n", maxDepth);


  }
}

TEST(LSGLRenderTest, TestBinaryPrefixTree60) {

  assert(sizeof(real) == sizeof(float));

  //int ns[] = {8, 16, 32, 1024, 1024*1024, 5555555, 123456789};
  //int ns[] = {8, 15, 31, 1023, 1000, 3211,  4000000, 12345678};
  int ns[] = {8, 15, 31, 1023, 1000, 3211};
  //int ns[] = {8, 16, 31};

  int niter = sizeof(ns) / sizeof(int);

  int seed = 123456;
  InitRandom(seed);

  for (int i = 0; i < niter; i++) {

    int npoints = ns[i];

    std::vector<float> points(3*npoints);
    std::vector<int> indices(npoints);
    std::vector<uint64_t> codes(npoints);

    for (int j = 0; j < npoints; j++) {
      points[3*j+0] = GenRandomReal();
      points[3*j+1] = GenRandomReal();
      points[3*j+2] = GenRandomReal();

      indices[j] = j;
    }

    real3 bmin(0.0, 0.0, 0.0);
    real3 bmax(1.0, 1.0, 1.0);

    CalculateMortonCodes60(&codes.at(0), &points.at(0), bmin, bmax, 0, npoints);

    std::vector<IndexKey60> keys(npoints);
    for (size_t i = 0; i < npoints; i++) {
      keys[i].index = indices[i];
      keys[i].code  = codes[i];
    }
    RadixSortByMortionCode60LSB(&keys.at(0), &keys.at(0) + npoints);
    //for (size_t i = 0; i < npoints; i++) {
    //  std::string k = BitString64(keys[i].code);
    //  printf("[%010d] = %010d(%s)\n", i, keys[i].code, k.c_str());
    //}

    // Tree size = npoints - 1
    std::vector<NodeInfo64> nodes(npoints-1);
    for (size_t i = 0; i < npoints - 1; i++) {
      nodes[i] = ConstructBinaryRadixTree60(&keys.at(0), i, npoints);
    }

    // Tree must visit each leaf node only once.
    {
      std::vector<int> counter(npoints);
      for (size_t i = 0; i < counter.size(); i++) {
        counter[i] = 0;
      }
    
      CountLeafs64(counter, nodes);

      //for (size_t i = 0; i < counter.size(); i++) {
      //    printf("count[%d] = %d\n", i, counter[i]);
      //}

      for (size_t i = 0; i < counter.size(); i++) {
        if (counter[i] != 1) {
          printf("i = %d, node = %d, leaf = %d, %d\n", i, nodes[i].index, nodes[i].leftType, nodes[i].rightType);
        }
        EXPECT_EQ(1, counter[i]);
      }
    }

    // Leaf tree must be accessible by traversing 
    {
      std::vector<int> counter(npoints);
      for (size_t i = 0; i < counter.size(); i++) {
        counter[i] = 0;
      }
    
      VisitAndCountLeafs64(counter, nodes, 0);

      for (size_t i = 0; i < counter.size(); i++) {
        if (counter[i] != 1) {
          printf("i = %d, node = %d, leaf = %d, %d\n", i, nodes[i].index, nodes[i].leftType, nodes[i].rightType);
        }
        EXPECT_EQ(1, counter[i]);
      }
    }

    // Check depth
    int maxDepth = 0;
    VisitAndCountDepth64(maxDepth, nodes, 0, 0);
    EXPECT_GT(60, maxDepth);
    //printf("maxDepth = %d\n", maxDepth);


  }
}

TEST(LSGLRenderTest, TestBinaryPrefixTreeSameKey30) {

  assert(sizeof(real) == sizeof(float));

  //int ns[] = {8, 16, 32, 1024, 1024*1024, 5555555, 123456789};
  //int ns[] = {8, 15, 31, 1023, 1000, 3211,  4000000, 12345678};
  int ns[] = {8, 15, 31, 1023, 1000, 3211};
  //int ns[] = {8, 16, 31};

  int niter = sizeof(ns) / sizeof(int);

  int seed = 123456;
  InitRandom(seed);

  for (int i = 0; i < niter; i++) {

    int npoints = ns[i];

    std::vector<float> points(3*npoints);
    std::vector<int> indices(npoints);
    std::vector<uint32_t> codes(npoints);

    for (int j = 0; j < npoints; j++) {
      points[3*j+0] = 0.123458;
      points[3*j+1] = 0.123458;
      points[3*j+2] = 0.123458;

      indices[j] = j;
    }

    real3 bmin(0.0, 0.0, 0.0);
    real3 bmax(1.0, 1.0, 1.0);

    CalculateMortonCodes30(&codes.at(0), &points.at(0), bmin, bmax, 0, npoints);

    std::vector<IndexKey30> keys(npoints);
    for (size_t i = 0; i < npoints; i++) {
      keys[i].index = indices[i];
      keys[i].code  = codes[i];
    }
    RadixSortByMortionCode30MSB(&keys.at(0), &keys.at(0) + npoints);

    // Tree size = npoints - 1
    std::vector<NodeInfo32> nodes(npoints-1);
    for (size_t i = 0; i < npoints - 1; i++) {
      nodes[i] = ConstructBinaryRadixTree30(&keys.at(0), i, npoints);
    }

    // Tree must visit each leaf node only once.
    {
      std::vector<int> counter(npoints);
      for (size_t i = 0; i < counter.size(); i++) {
        counter[i] = 0;
      }
    
      CountLeafs32(counter, nodes);

      //for (size_t i = 0; i < counter.size(); i++) {
      //    printf("count[%d] = %d\n", i, counter[i]);
      //}

      for (size_t i = 0; i < counter.size(); i++) {
        if (counter[i] != 1) {
          printf("i = %d, node = %d, leaf = %d, %d\n", i, nodes[i].childIndex, nodes[i].leftType, nodes[i].rightType);
        }
        EXPECT_EQ(1, counter[i]);
      }
    }

    // Leaf tree must be accessible by traversing 
    {
      std::vector<int> counter(npoints);
      for (size_t i = 0; i < counter.size(); i++) {
        counter[i] = 0;
      }
    
      VisitAndCountLeafs32(counter, nodes, 0);

      for (size_t i = 0; i < counter.size(); i++) {
        if (counter[i] != 1) {
          printf("i = %d, node = %d, leaf = %d, %d\n", i, nodes[i].childIndex, nodes[i].leftType, nodes[i].rightType);
        }
        EXPECT_EQ(1, counter[i]);
      }
    }

    // Check depth
    int maxDepth = 0;
    VisitAndCountDepth32(maxDepth, nodes, 0, 0);
    EXPECT_GT(30, maxDepth);
    //printf("maxDepth = %d\n", maxDepth);


  }
}

TEST(LSGLRenderTest, TestBinaryPrefixTreeSameKey60) {

  assert(sizeof(real) == sizeof(float));

  //int ns[] = {8, 16, 32, 1024, 1024*1024, 5555555, 123456789};
  //int ns[] = {8, 15, 31, 1023, 1000, 3211,  4000000, 12345678};
  int ns[] = {8, 15, 31, 1023, 1000, 3211};
  //int ns[] = {8, 16, 31};

  int niter = sizeof(ns) / sizeof(int);

  int seed = 123456;
  InitRandom(seed);

  for (int i = 0; i < niter; i++) {

    int npoints = ns[i];

    std::vector<float> points(3*npoints);
    std::vector<int> indices(npoints);
    std::vector<uint64_t> codes(npoints);

    for (int j = 0; j < npoints; j++) {
      points[3*j+0] = 0.123458;
      points[3*j+1] = 0.123458;
      points[3*j+2] = 0.123458;

      indices[j] = j;
    }

    real3 bmin(0.0, 0.0, 0.0);
    real3 bmax(1.0, 1.0, 1.0);

    CalculateMortonCodes60(&codes.at(0), &points.at(0), bmin, bmax, 0, npoints);

    std::vector<IndexKey60> keys(npoints);
    for (size_t i = 0; i < npoints; i++) {
      keys[i].index = indices[i];
      keys[i].code  = codes[i];
    }
    RadixSortByMortionCode60MSB(&keys.at(0), &keys.at(0) + npoints);

    // Tree size = npoints - 1
    std::vector<NodeInfo64> nodes(npoints-1);
    for (size_t i = 0; i < npoints - 1; i++) {
      nodes[i] = ConstructBinaryRadixTree60(&keys.at(0), i, npoints);
    }

    // Tree must visit each leaf node only once.
    {
      std::vector<int> counter(npoints);
      for (size_t i = 0; i < counter.size(); i++) {
        counter[i] = 0;
      }
    
      CountLeafs64(counter, nodes);

      for (size_t i = 0; i < counter.size(); i++) {
        if (counter[i] != 1) {
          printf("i = %d, node = %d, leaf = %d, %d\n", i, nodes[i].index, nodes[i].leftType, nodes[i].rightType);
        }
        EXPECT_EQ(1, counter[i]);
      }
    }

    // Leaf tree must be accessible by traversing 
    {
      std::vector<int> counter(npoints);
      for (size_t i = 0; i < counter.size(); i++) {
        counter[i] = 0;
      }
    
      VisitAndCountLeafs64(counter, nodes, 0);

      for (size_t i = 0; i < counter.size(); i++) {
        if (counter[i] != 1) {
          printf("i = %d, node = %d, leaf = %d, %d\n", i, nodes[i].index, nodes[i].leftType, nodes[i].rightType);
        }
        EXPECT_EQ(1, counter[i]);
      }
    }

    // Check depth
    int maxDepth = 0;
    VisitAndCountDepth64(maxDepth, nodes, 0, 0);
    EXPECT_GT(60, maxDepth);
    //printf("maxDepth = %d\n", maxDepth);


  }
}
