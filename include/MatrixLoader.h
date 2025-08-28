#ifndef MATRIX_LOADER_H
#define MATRIX_LOADER_H

#include <iostream>
#include <vector>
#include <string>
#include <variant>

namespace MatrixProject {

using Element = std::variant<double,std::vector<double>>;
using Matrix = std::vector<std::vector<Element>>;
using MatrixGroup = std::vector<Matrix>;
using MatrixGroups = std::vector<MatrixGroup>;

/**
 * @brief 读取矩阵数据从文本文件
 * @param filename 文件名
 * @return 返回矩阵组的向量
 */
MatrixGroups read_matrix_from_txt(const std::string& filename);

/**
 * @brief 对矩阵元素进行排序
 * @param matrix_groups 矩阵组
 * @return 排序后的矩阵组
 */
MatrixGroups sort_matrix_elements(const MatrixGroups& matrix_groups);

/**
 * @brief 对矩阵组中的矩阵进行去重
 * @param matrix_groups 矩阵组
 * @return 去重后的矩阵组
 */
MatrixGroups remove_duplicate_matrices(const MatrixGroups& matrix_groups);

} // namespace MatrixProject

#endif // MATRIX_LOADER_H