// [leukoquant patch] See _reg_maths_eigen_compat.h for why this file exists.
#include "_reg_maths_eigen_compat.h"

#include <cstdio>
#include <cstdlib>
#include <algorithm>

// Eigen headers are in there because of the nvcc preprocessing step
// [leukoquant patch] The current vendored NiftyReg's third-party/Eigen/ is a
// flattened distribution -- MatrixFunctions (an "unsupported" module in a
// standard Eigen release) lives directly under Eigen/, not
// unsupported/Eigen/, so this is Eigen/MatrixFunctions rather than the
// unsupported/Eigen/MatrixFunctions path the old code used.
#include "Eigen/Core"
#include "Eigen/SVD"
#include "Eigen/MatrixFunctions"

/* *************************************************************** */
// Self-contained replacements for NiftyReg helpers this file used to rely
// on (reg_print_fct_error/reg_print_msg_error/reg_exit, reg_mat33_to_nan,
// reg_mat44_add, reg_mat44_mul) that no longer exist in the current
// vendored NiftyReg -- kept local so this file has no further dependency
// on symbols that could disappear in a future NiftyReg refactor.
static void gif_maths_fatal_error(const char *function, const char *message) {
   fprintf(stderr, "[%s] ERROR: %s\n", function, message);
   exit(1);
}
static void gif_mat33_to_nan(mat33 *m) {
   for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
         m->m[i][j] = std::numeric_limits<float>::quiet_NaN();
}
/* *************************************************************** */
void reg_print_msg_error(const char *msg) {
   fprintf(stderr, "[NiftyReg ERROR] %s\n", msg);
}
/* *************************************************************** */
mat44 reg_mat44_mul(mat44 const* A, mat44 const* B) {
   mat44 R;
   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         R.m[i][j] = static_cast<float>(
            static_cast<double>(A->m[i][0]) * static_cast<double>(B->m[0][j]) +
            static_cast<double>(A->m[i][1]) * static_cast<double>(B->m[1][j]) +
            static_cast<double>(A->m[i][2]) * static_cast<double>(B->m[2][j]) +
            static_cast<double>(A->m[i][3]) * static_cast<double>(B->m[3][j]));
   return R;
}
/* *************************************************************** */
mat44 reg_mat44_mul(mat44 const* A, double scalar) {
   mat44 out;
   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         out.m[i][j] = static_cast<float>(A->m[i][j] * scalar);
   return out;
}
/* *************************************************************** */
mat44 reg_mat44_add(mat44 const* A, mat44 const* B) {
   mat44 R;
   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         R.m[i][j] = static_cast<float>(static_cast<double>(A->m[i][j]) + static_cast<double>(B->m[i][j]));
   return R;
}
/* *************************************************************** */
mat44 reg_mat44_minus(mat44 const* A, mat44 const* B) {
   mat44 R;
   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
         R.m[i][j] = static_cast<float>(static_cast<double>(A->m[i][j]) - static_cast<double>(B->m[i][j]));
   return R;
}
/* *************************************************************** */
void reg_mat44_disp(mat44 *mat, char *title) {
   printf("%s:\n%.7g\t%.7g\t%.7g\t%.7g\n%.7g\t%.7g\t%.7g\t%.7g\n%.7g\t%.7g\t%.7g\t%.7g\n%.7g\t%.7g\t%.7g\t%.7g\n", title,
          mat->m[0][0], mat->m[0][1], mat->m[0][2], mat->m[0][3],
          mat->m[1][0], mat->m[1][1], mat->m[1][2], mat->m[1][3],
          mat->m[2][0], mat->m[2][1], mat->m[2][2], mat->m[2][3],
          mat->m[3][0], mat->m[3][1], mat->m[3][2], mat->m[3][3]);
}
/* *************************************************************** */
/** @brief SVD
* @param in input matrix to decompose - in place
* @param size_m row
* @param size_n colomn
* @param w diagonal term
* @param v rotation part
*/
template<class T>
void svd(T **in, size_t size_m, size_t size_n, T * w, T **v) {
   if (size_m == 0 || size_n == 0)
      gif_maths_fatal_error("svd", "The specified matrix is empty");

   size_t sm, sn, sn2;
   Eigen::MatrixXd m(size_m, size_n);

   for (sm = 0; sm < size_m; sm++)
      for (sn = 0; sn < size_n; sn++)
         m(sm, sn) = static_cast<double>(in[sm][sn]);

   Eigen::JacobiSVD<Eigen::MatrixXd> svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);

   for (sn = 0; sn < size_n; sn++) {
      w[sn] = static_cast<T>(svd.singularValues()(sn));
      for (sn2 = 0; sn2 < size_n; sn2++)
         v[sn2][sn] = static_cast<T>(svd.matrixV()(sn2, sn));
      for (sm = 0; sm < size_m; sm++)
         in[sm][sn] = static_cast<T>(svd.matrixU()(sm, sn));
   }
}
template void svd<float>(float **in, size_t m, size_t n, float * w, float **v);
template void svd<double>(double **in, size_t m, size_t n, double * w, double **v);
/* *************************************************************** */
/**
* @brief SVD
* @param in input matrix to decompose
* @param size_m row
* @param size_n colomn
* @param U unitary matrices
* @param S diagonal matrix
* @param V unitary matrices
*  X = U*S*V'
*/
template<class T>
void svd(T **in, size_t size_m, size_t size_n, T ***U, T ***S, T ***V) {
   if (in == nullptr)
      gif_maths_fatal_error("svd", "The specified matrix is empty");

   size_t sm, sn, min_dim, i, j;
   Eigen::MatrixXd m(size_m, size_n);

   for (sm = 0; sm < size_m; sm++)
      for (sn = 0; sn < size_n; sn++)
         m(sm, sn) = static_cast<double>(in[sm][sn]);

   Eigen::JacobiSVD<Eigen::MatrixXd> svd(m, Eigen::ComputeThinU | Eigen::ComputeThinV);

   min_dim = std::min(size_m, size_n);
   for (i = 0; i < min_dim; i++) {
      for (j = 0; j < min_dim; j++) {
         if (i == j)
            (*S)[i][j] = static_cast<T>(svd.singularValues()(i));
         else
            (*S)[i][j] = 0;
      }
   }

   if (size_m > size_n) {
      for (i = 0; i < min_dim; i++)
         for (j = 0; j < min_dim; j++)
            (*V)[i][j] = static_cast<T>(svd.matrixV()(i, j));
      for (i = 0; i < size_m; i++)
         for (j = 0; j < size_n; j++)
            (*U)[i][j] = static_cast<T>(svd.matrixU()(i, j));
   } else {
      for (i = 0; i < min_dim; i++)
         for (j = 0; j < min_dim; j++)
            (*U)[i][j] = static_cast<T>(svd.matrixU()(i, j));
      for (i = 0; i < size_n; i++)
         for (j = 0; j < size_m; j++)
            (*V)[i][j] = static_cast<T>(svd.matrixV()(i, j));
   }
}
template void svd<float>(float **in, size_t size_m, size_t size_n, float ***U, float ***S, float ***V);
template void svd<double>(double **in, size_t size_m, size_t size_n, double ***U, double ***S, double ***V);
/* *************************************************************** */
template<class T>
T reg_matrix2DDet(T** mat, size_t m, size_t n) {
   if (m != n) {
      char text[255];
      sprintf(text, "The matrix have to be square: [%zu %zu]", m, n);
      gif_maths_fatal_error("reg_matrix2DDet", text);
   }
   double res;
   if (m == 2) {
      res = static_cast<double>(mat[0][0]) * static_cast<double>(mat[1][1]) - static_cast<double>(mat[1][0]) * static_cast<double>(mat[0][1]);
   } else if (m == 3) {
      res = (static_cast<double>(mat[0][0]) * (static_cast<double>(mat[1][1]) * static_cast<double>(mat[2][2]) - static_cast<double>(mat[1][2]) * static_cast<double>(mat[2][1]))) -
            (static_cast<double>(mat[0][1]) * (static_cast<double>(mat[1][0]) * static_cast<double>(mat[2][2]) - static_cast<double>(mat[1][2]) * static_cast<double>(mat[2][0]))) +
            (static_cast<double>(mat[0][2]) * (static_cast<double>(mat[1][0]) * static_cast<double>(mat[2][1]) - static_cast<double>(mat[1][1]) * static_cast<double>(mat[2][0])));
   } else {
      Eigen::MatrixXd eigenRes(m, n);
      for (size_t i = 0; i < m; i++)
         for (size_t j = 0; j < n; j++)
            eigenRes(i, j) = static_cast<double>(mat[i][j]);
      res = eigenRes.determinant();
   }
   return static_cast<T>(res);
}
template float reg_matrix2DDet<float>(float** mat, size_t m, size_t n);
template double reg_matrix2DDet<double>(double** mat, size_t m, size_t n);
/* *************************************************************** */
mat44 reg_mat44_sqrt(mat44 const* mat) {
   mat44 X;
   Eigen::Matrix4d m;
   for (size_t i = 0; i < 4; ++i)
      for (size_t j = 0; j < 4; ++j)
         m(i, j) = static_cast<double>(mat->m[i][j]);
   m = m.sqrt();
   for (size_t i = 0; i < 4; ++i)
      for (size_t j = 0; j < 4; ++j)
         X.m[i][j] = static_cast<float>(m(i, j));
   return X;
}
/* *************************************************************** */
void reg_mat33_expm(mat33 *in_tensor) {
   int sm, sn;
   Eigen::Matrix3d tensor;
   for (sm = 0; sm < 3; sm++) {
      for (sn = 0; sn < 3; sn++) {
         float val = in_tensor->m[sm][sn];
         if (val != val) return;
         tensor(sm, sn) = static_cast<double>(val);
      }
   }
   tensor = tensor.exp();
   for (sm = 0; sm < 3; sm++)
      for (sn = 0; sn < 3; sn++)
         in_tensor->m[sm][sn] = static_cast<float>(tensor(sm, sn));
}
/* *************************************************************** */
mat44 reg_mat44_expm(mat44 const* mat) {
   mat44 X;
   Eigen::Matrix4d m;
   for (size_t i = 0; i < 4; ++i)
      for (size_t j = 0; j < 4; ++j)
         m(i, j) = static_cast<double>(mat->m[i][j]);
   m = m.exp();
   for (size_t i = 0; i < 4; ++i)
      for (size_t j = 0; j < 4; ++j)
         X.m[i][j] = static_cast<float>(m(i, j));
   return X;
}
/* *************************************************************** */
void reg_mat33_logm(mat33 *in_tensor) {
   int sm, sn;
   Eigen::Matrix3d tensor;
   bool all_zeros = true;
   double det = 0;
   for (sm = 0; sm < 3; sm++) {
      for (sn = 0; sn < 3; sn++) {
         float val = in_tensor->m[sm][sn];
         if (val != 0.f) all_zeros = false;
         if (val != val) return;
         tensor(sm, sn) = static_cast<double>(val);
      }
   }
   det = tensor.determinant();
   if (all_zeros || det == 0) {
      gif_mat33_to_nan(in_tensor);
      return;
   }
   tensor = tensor.log();
   for (sm = 0; sm < 3; sm++)
      for (sn = 0; sn < 3; sn++)
         in_tensor->m[sm][sn] = static_cast<float>(tensor(sm, sn));
}
/* *************************************************************** */
mat44 reg_mat44_logm(mat44 const* mat) {
   mat44 X;
   Eigen::Matrix4d m;
   for (size_t i = 0; i < 4; ++i)
      for (size_t j = 0; j < 4; ++j)
         m(i, j) = static_cast<double>(mat->m[i][j]);
   m = m.log();
   for (size_t i = 0; i < 4; ++i)
      for (size_t j = 0; j < 4; ++j)
         X.m[i][j] = static_cast<float>(m(i, j));
   return X;
}
/* *************************************************************** */
mat44 reg_mat44_inv(mat44 const* mat) {
   mat44 out;
   Eigen::Matrix4d m, m_inv;
   for (size_t i = 0; i < 4; ++i)
      for (size_t j = 0; j < 4; ++j)
         m(i, j) = static_cast<double>(mat->m[i][j]);
   m_inv = m.inverse();
   for (size_t i = 0; i < 4; ++i)
      for (size_t j = 0; j < 4; ++j)
         out.m[i][j] = static_cast<float>(m_inv(i, j));
   return out;
}
/* *************************************************************** */
mat44 reg_mat44_avg2(mat44 const* A, mat44 const* B) {
   mat44 logA = reg_mat44_logm(A);
   mat44 logB = reg_mat44_logm(B);
   for (int i = 0; i < 4; ++i) {
      logA.m[3][i] = 0.f;
      logB.m[3][i] = 0.f;
   }
   mat44 sum = reg_mat44_add(&logA, &logB);
   mat44 avg = reg_mat44_mul(&sum, 0.5);
   return reg_mat44_expm(&avg);
}
