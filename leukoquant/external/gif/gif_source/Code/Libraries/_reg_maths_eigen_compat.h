#pragma once

#include "nifti1_io.h"

/* *************************************************************** */
/* [leukoquant patch] The version of NiftyReg vendored under
 * leukoquant/external/niftyreg/source/ (post-#130 "reg_resample cuda
 * enabling" refactor) removed _reg_maths_eigen.h/.cpp entirely -- these
 * Eigen-backed matrix functions (SVD, matrix exp/log/inverse/sqrt) no
 * longer exist anywhere in NiftyReg. GIF's own source (_seg_GIF.cpp) still
 * calls reg_mat44_expm(), so this file ports the full old API locally,
 * self-contained (not depending on reg_print_fct_error/reg_exit/
 * reg_mat33_to_nan/reg_mat44_add/reg_mat44_mul, which were also removed).
 * Ported from the pre-refactor NiftyReg source (reg-lib/cpu/
 * _reg_maths_eigen.h/.cpp), functionally unchanged.
 * *************************************************************** */

extern "C++" template <class T>
void svd(T **in, size_t m, size_t n, T *w, T **v);
/* *************************************************************** */
extern "C++" template <class T>
void svd(T **in, size_t m, size_t n, T ***U, T ***S, T ***V);
/* *************************************************************** */
extern "C++" template<class T>
T reg_matrix2DDet(T** mat, size_t m, size_t n);
/* *************************************************************** */
/** @brief Compute the inverse of a  4-by-4 matrix
*/
mat44 reg_mat44_inv(mat44 const* mat);
/* *************************************************************** */
/** @brief Compute the square root of a 4-by-4 matrix
*/
mat44 reg_mat44_sqrt(mat44 const* mat);
/* *************************************************************** */
/** @brief Compute the exp of a 3-by-3 matrix
*/
void reg_mat33_expm(mat33 *in_tensor);
/* *************************************************************** */
/** @brief Compute the exp of a 4-by-4 matrix
*/
mat44 reg_mat44_expm(const mat44 *mat);
/* *************************************************************** */
/** @brief Compute the log of a 3-by-3 matrix
*/
void reg_mat33_logm(mat33 *in_tensor);
/* *************************************************************** */
/** @brief Compute the log of a 4-by-4 matrix
*/
mat44 reg_mat44_logm(const mat44 *mat);
/* *************************************************************** */
/** @brief Compute the average of two matrices using a log-euclidean
* framework
*/
mat44 reg_mat44_avg2(mat44 const* A, mat44 const* b);
/* *************************************************************** */
/* Below: plain (non-Eigen) mat44 arithmetic/display helpers that also no
 * longer exist in the current vendored NiftyReg (were in its old
 * _reg_maths.h, not _reg_maths_eigen.h, but removed the same way). Grouped
 * here rather than a second compat file since GIF needs both together. */
mat44 reg_mat44_mul(mat44 const* A, mat44 const* B);
mat44 reg_mat44_mul(mat44 const* A, double scalar);
mat44 reg_mat44_add(mat44 const* A, mat44 const* B);
mat44 reg_mat44_minus(mat44 const* A, mat44 const* B);
void reg_mat44_disp(mat44 *mat, char *title);
/* *************************************************************** */
/* Two more small NiftyReg utility functions GIF calls that also no longer
 * exist in the current vendored NiftyReg. reg_pow2 is inlined here (not in
 * the .cpp) since it's a template -- needs to be visible at each call site's
 * instantiation, not just explicitly instantiated for one type. */
void reg_print_msg_error(const char *msg);
template<class T> inline T reg_pow2(T a) { return a * a; }
