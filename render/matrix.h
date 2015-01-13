#ifndef __RENDER_MATRIX_H__
#define __RENDER_MATRIX_H__

#include "vector3.h"

namespace lsgl {
namespace render {

// 4x4 matrix
template <typename T> class Matrix {
public:
  Matrix();
  ~Matrix();

  static void Print(T m[4][4]) {
    for (int i = 0; i < 4; i++) {
      printf("m[%d] = %f, %f, %f, %f\n", i, m[i][0], m[i][1], m[i][2], m[i][3]);
    }
  }

  static void Identity(T m[4][4]) {
    m[0][0] = 1.0;
    m[0][1] = 0.0;
    m[0][2] = 0.0;
    m[0][3] = 0.0;
    m[1][0] = 0.0;
    m[1][1] = 1.0;
    m[1][2] = 0.0;
    m[1][3] = 0.0;
    m[2][0] = 0.0;
    m[2][1] = 0.0;
    m[2][2] = 1.0;
    m[2][3] = 0.0;
    m[3][0] = 0.0;
    m[3][1] = 0.0;
    m[3][2] = 0.0;
    m[3][3] = 1.0;
  }

  static void Inverse(T m[4][4]) {
    /*
     * codes from intel web
     * cramer's rule version
     */
    int i, j;
    T tmp[12];  /* tmp array for pairs */
    T tsrc[16]; /* array of transpose source matrix */
    T det;      /* determinant */

    /* transpose matrix */
    for (i = 0; i < 4; i++) {
      tsrc[i] = m[i][0];
      tsrc[i + 4] = m[i][1];
      tsrc[i + 8] = m[i][2];
      tsrc[i + 12] = m[i][3];
    }

    /* calculate pair for first 8 elements(cofactors) */
    tmp[0] = tsrc[10] * tsrc[15];
    tmp[1] = tsrc[11] * tsrc[14];
    tmp[2] = tsrc[9] * tsrc[15];
    tmp[3] = tsrc[11] * tsrc[13];
    tmp[4] = tsrc[9] * tsrc[14];
    tmp[5] = tsrc[10] * tsrc[13];
    tmp[6] = tsrc[8] * tsrc[15];
    tmp[7] = tsrc[11] * tsrc[12];
    tmp[8] = tsrc[8] * tsrc[14];
    tmp[9] = tsrc[10] * tsrc[12];
    tmp[10] = tsrc[8] * tsrc[13];
    tmp[11] = tsrc[9] * tsrc[12];

    /* calculate first 8 elements(cofactors) */
    m[0][0] = tmp[0] * tsrc[5] + tmp[3] * tsrc[6] + tmp[4] * tsrc[7];
    m[0][0] -= tmp[1] * tsrc[5] + tmp[2] * tsrc[6] + tmp[5] * tsrc[7];
    m[0][1] = tmp[1] * tsrc[4] + tmp[6] * tsrc[6] + tmp[9] * tsrc[7];
    m[0][1] -= tmp[0] * tsrc[4] + tmp[7] * tsrc[6] + tmp[8] * tsrc[7];
    m[0][2] = tmp[2] * tsrc[4] + tmp[7] * tsrc[5] + tmp[10] * tsrc[7];
    m[0][2] -= tmp[3] * tsrc[4] + tmp[6] * tsrc[5] + tmp[11] * tsrc[7];
    m[0][3] = tmp[5] * tsrc[4] + tmp[8] * tsrc[5] + tmp[11] * tsrc[6];
    m[0][3] -= tmp[4] * tsrc[4] + tmp[9] * tsrc[5] + tmp[10] * tsrc[6];
    m[1][0] = tmp[1] * tsrc[1] + tmp[2] * tsrc[2] + tmp[5] * tsrc[3];
    m[1][0] -= tmp[0] * tsrc[1] + tmp[3] * tsrc[2] + tmp[4] * tsrc[3];
    m[1][1] = tmp[0] * tsrc[0] + tmp[7] * tsrc[2] + tmp[8] * tsrc[3];
    m[1][1] -= tmp[1] * tsrc[0] + tmp[6] * tsrc[2] + tmp[9] * tsrc[3];
    m[1][2] = tmp[3] * tsrc[0] + tmp[6] * tsrc[1] + tmp[11] * tsrc[3];
    m[1][2] -= tmp[2] * tsrc[0] + tmp[7] * tsrc[1] + tmp[10] * tsrc[3];
    m[1][3] = tmp[4] * tsrc[0] + tmp[9] * tsrc[1] + tmp[10] * tsrc[2];
    m[1][3] -= tmp[5] * tsrc[0] + tmp[8] * tsrc[1] + tmp[11] * tsrc[2];

    /* calculate pairs for second 8 elements(cofactors) */
    tmp[0] = tsrc[2] * tsrc[7];
    tmp[1] = tsrc[3] * tsrc[6];
    tmp[2] = tsrc[1] * tsrc[7];
    tmp[3] = tsrc[3] * tsrc[5];
    tmp[4] = tsrc[1] * tsrc[6];
    tmp[5] = tsrc[2] * tsrc[5];
    tmp[6] = tsrc[0] * tsrc[7];
    tmp[7] = tsrc[3] * tsrc[4];
    tmp[8] = tsrc[0] * tsrc[6];
    tmp[9] = tsrc[2] * tsrc[4];
    tmp[10] = tsrc[0] * tsrc[5];
    tmp[11] = tsrc[1] * tsrc[4];

    /* calculate second 8 elements(cofactors) */
    m[2][0] = tmp[0] * tsrc[13] + tmp[3] * tsrc[14] + tmp[4] * tsrc[15];
    m[2][0] -= tmp[1] * tsrc[13] + tmp[2] * tsrc[14] + tmp[5] * tsrc[15];
    m[2][1] = tmp[1] * tsrc[12] + tmp[6] * tsrc[14] + tmp[9] * tsrc[15];
    m[2][1] -= tmp[0] * tsrc[12] + tmp[7] * tsrc[14] + tmp[8] * tsrc[15];
    m[2][2] = tmp[2] * tsrc[12] + tmp[7] * tsrc[13] + tmp[10] * tsrc[15];
    m[2][2] -= tmp[3] * tsrc[12] + tmp[6] * tsrc[13] + tmp[11] * tsrc[15];
    m[2][3] = tmp[5] * tsrc[12] + tmp[8] * tsrc[13] + tmp[11] * tsrc[14];
    m[2][3] -= tmp[4] * tsrc[12] + tmp[9] * tsrc[13] + tmp[10] * tsrc[14];
    m[3][0] = tmp[2] * tsrc[10] + tmp[5] * tsrc[11] + tmp[1] * tsrc[9];
    m[3][0] -= tmp[4] * tsrc[11] + tmp[0] * tsrc[9] + tmp[3] * tsrc[10];
    m[3][1] = tmp[8] * tsrc[11] + tmp[0] * tsrc[8] + tmp[7] * tsrc[10];
    m[3][1] -= tmp[6] * tsrc[10] + tmp[9] * tsrc[11] + tmp[1] * tsrc[8];
    m[3][2] = tmp[6] * tsrc[9] + tmp[11] * tsrc[11] + tmp[3] * tsrc[8];
    m[3][2] -= tmp[10] * tsrc[11] + tmp[2] * tsrc[8] + tmp[7] * tsrc[9];
    m[3][3] = tmp[10] * tsrc[10] + tmp[4] * tsrc[8] + tmp[9] * tsrc[9];
    m[3][3] -= tmp[8] * tsrc[9] + tmp[11] * tsrc[0] + tmp[5] * tsrc[8];

    /* calculate determinant */
    det = tsrc[0] * m[0][0] + tsrc[1] * m[0][1] + tsrc[2] * m[0][2] +
          tsrc[3] * m[0][3];

    /* calculate matrix inverse */
    det = 1.0 / det;

    for (j = 0; j < 4; j++) {
      for (i = 0; i < 4; i++) {
        m[j][i] *= det;
      }
    }
  }

  static void Transpose(T m[4][4]) {
    T t[4][4];

    // Transpose
    for (int j = 0; j < 4; j++) {
      for (int i = 0; i < 4; i++) {
        t[j][i] = m[i][j];
      }
    }

    // Copy
    for (int j = 0; j < 4; j++) {
      for (int i = 0; i < 4; i++) {
        m[j][i] = t[j][i];
      }
    }
  }

  static void Mult(T dst[4][4], const T m0[4][4], const T m1[4][4]) {
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
        dst[i][j] = 0;
        for (int k = 0; k < 4; ++k) {
          dst[i][j] += m0[k][j] * m1[i][k];
        }
      }
    }
  }

  static void MultV(T dst[3], const T m[4][4], const T v[3]) {
    T tmp[3];
    tmp[0] = m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2] + m[3][0];
    tmp[1] = m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2] + m[3][1];
    tmp[2] = m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2] + m[3][2];
    dst[0] = tmp[0];
    dst[1] = tmp[1];
    dst[2] = tmp[2];
  }

  static void MultV(vector3 &dst, const T m[4][4], const vector3 &v) {
    T tmp[3];
    tmp[0] = m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2] + m[3][0];
    tmp[1] = m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2] + m[3][1];
    tmp[2] = m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2] + m[3][2];
    dst[0] = tmp[0];
    dst[1] = tmp[1];
    dst[2] = tmp[2];
  }
};

typedef Matrix<float> Matrixf;
typedef Matrix<double> Matrixd;

} // namespace render
} // namespace lsgl

#endif // __RENDER_MATRIX_H__
