#ifndef MATRIX_AB_H
#define MATRIX_AB_H

#include "MatrixLoader.h"
#include "MatrixTools.h"
#include <vector>

namespace MatrixProject {

// 辅助函数：检查 Element 是否为零
bool is_zero(const Element& elem);

/**
 * @brief 检查给定元素是否存在于矩阵的某一列中，包括列表中的元素。
 * @param matrix 输入矩阵
 * @param element 要检查的元素
 * @param col 列索引
 * @return 如果存在，返回 true；否则返回 false
 */
bool element_in_column(const Matrix& matrix, double element, size_t col);

/**
 * @brief 检查元素是否不在给定的矩阵位置中，包括列表中的元素。
 * @param position 矩阵的特定位置元素
 * @param element 要检查的元素
 * @return 如果元素不在位置中，返回 true；否则返回 false
 */
bool element_not_in_position(const Element& position, double element);

/**
 * @brief 合并多个矩阵，生成矩阵 A 和 B
 * @param matrices 输入的矩阵列表
 * @param needPrint 是否需要打印结果
 * @return 一个 pair，包含矩阵 A 和矩阵 B
 */
// std::pair<Matrix, Matrix> merge_matrices_custom(const std::vector<Matrix>& matrices, bool needPrint);
std::tuple<Matrix, Matrix> merge_matrices_custom(const std::vector<Matrix>& matrices, bool needPrint);

} // namespace MatrixProject

#endif // MATRIX_AB_H