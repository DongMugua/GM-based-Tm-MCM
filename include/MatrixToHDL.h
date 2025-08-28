//
// Created by 孙皓 on 2025/3/14.
//

#ifndef MATRIXTOHDL_H
#define MATRIXTOHDL_H

#include <string>
#include <vector>
#include <variant>
#include "MatrixLoader.h"

namespace MatrixProject {

/**
 * @brief 从指定的 txt 文件中解析出：
 *        1) Selected Matrices（可能包含多个矩阵）
 *        2) 矩阵 A
 *        3) 矩阵 B
 *
 * @param filename           文件路径
 * @param selected_matrices  输出的 Selected Matrices
 * @param matrixA            输出的 Matrix A
 * @param matrixB            输出的 Matrix B
 * @return true  解析成功
 * @return false 解析失败（无法打开文件等）
 */
bool parse_matrices_from_file(const std::string& filename,
                              std::vector<Matrix>& selected_matrices,
                              Matrix& matrixA,
                              Matrix& matrixB);
}

#endif // MATRIXTOHDL_H