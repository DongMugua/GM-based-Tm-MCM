#ifndef MATRIX_COST_H
#define MATRIX_COST_H

#include "MatrixLoader.h"
#include <vector>

namespace MatrixProject {

/**
 * @brief 计算矩阵非对角线上的非零元素数量（考虑重复元素）
 * @param matrix 输入矩阵
 * @return 非零元素数量的向量
 */
std::vector<int> count_non_zero_elements_off_diagonal(const Matrix& matrix);

/**
 * @brief 获取元素的绝对值最大值或原值中的绝对值最大值
 * @param element 输入元素（可能是数字或向量）
 * @param return_abs 是否返回绝对值，默认为 true
 * @return 元素的绝对值最大值或原值
 */
double get_max_value(const Element& element, bool return_abs = true);
double get_min_value(const Element& element, bool return_abs = true);

/**
 * @brief 对矩阵进行变换操作
 * @param matrix 输入矩阵
 * @return 变换后的矩阵
 */
Matrix transform_matrix(const Matrix& matrix);

/**
 * @brief 获取矩阵每列的最大值向量
 * @param matrix 输入矩阵
 * @return 每列最大值的向量
 */
std::vector<double> get_column_vector(const Matrix& matrix);

/**
 * @brief 计算 MUX 的成本和数量
 * @param matrixA 矩阵 A
 * @param matrixB 矩阵 B
 * @param mac_vector MAC 向量
 * @return 一个 pair，包含总成本和 MUX 数量
 */
std::pair<double, double> cost_matrices_mux(const Matrix& matrixA, const Matrix& matrixB, const std::vector<double>& mac_vector);

/**
 * @brief 分析加法器的成本
 * @param matrixA 矩阵 A
 * @param matrixB 矩阵 B
 * @return 加法器成本的向量
 */
std::vector<int> cost_matrices_adders(const Matrix& matrixA, const Matrix& matrixB);

/**
 * @brief 计算矩阵的总成本和 MUX 数量
 * @param matrixA 矩阵 A
 * @param matrixB 矩阵 B
 * @return 一个 pair，包含总成本和 MUX 数量
 */
std::pair<double, double> cost_analysis(const Matrix& matrixA, const Matrix& matrixB);

} // namespace MatrixProject

#endif // MATRIX_COST_H