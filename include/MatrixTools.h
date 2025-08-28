#ifndef MATRIX_TOOLS_H
#define MATRIX_TOOLS_H

#include "MatrixLoader.h"

#include <iostream>

namespace MatrixProject {
size_t element_width(const std::variant<double, std::vector<double>>& elem);

/**
 * @brief 打印矩阵
 * @param matrix 矩阵
 */
void print_matrix(const Matrix& matrix);
void print_matrix(const MatrixGroup& MatrixGroup);
void print_matrix(const MatrixGroups& MatrixGroups);
void print_matrix(const Matrix& matrix, std::ostream& out);

/**
 * @brief 移除矩阵中的重复向量元素
 * @param matrix 矩阵
 * @return 处理后的矩阵
 */
Matrix remove_duplicates_from_vectors(const Matrix& matrix);

/**
 * @brief 确保矩阵中的所有元素都是列表格式
 * @param matrix 矩阵
 * @return 处理后的矩阵
 */
Matrix ensure_list_format(const Matrix& matrix);

/**
 * @brief 将矩阵组中的所有矩阵填充到相同大小
 * @param matrix_groups 矩阵组
 * @return 填充后的矩阵组
 */
MatrixGroups pad_all_groups(const MatrixGroups& matrix_groups);

/**
 * @brief 获取所有矩阵的最大方阵尺寸
 * @param matrix_groups 矩阵组
 * @return 最大方阵尺寸
 */
int get_max_square_matrix_size(const MatrixGroups& matrix_groups);

/**
 * @brief 将矩阵填充为指定大小的方阵
 * @param matrix 矩阵
 * @param size 目标尺寸
 * @return 填充后的矩阵
 */
Matrix pad_square_matrix(const Matrix& matrix, int size);

/**
 * @brief 合并两个矩阵
 * @param matrixA 矩阵 A
 * @param matrixB 矩阵 B
 * @return 合并后的矩阵
 */
Matrix merge_matrices(const Matrix& matrixA, const Matrix& matrixB);

/**
 * @brief 合并两个矩阵
 * @param matrix
 * @param Size
 * @return 合并后的矩阵
 */
std::vector<Matrix> iterative_pad_matrix(const Matrix& matrix, int max_size);

// 检查两个 double 是否相等（考虑容差）
bool is_close(double a, double b, double tol);

// 检查两个 Element 是否相等
bool elements_equal(const Element& e1, const Element& e2, double tol);

// 检查矩阵是否存在于矩阵组中
bool matrix_exists(const Matrix& matrix, const MatrixGroup& matrix_group, double tol);

} // namespace MatrixProject

#endif // MATRIX_TOOLS_H