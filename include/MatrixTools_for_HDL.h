#ifndef MATRIX_TOOLS_FOR_HDL_H
#define MATRIX_TOOLS_FOR_HDL_H

#include "MatrixLoader.h"
#include <vector>
#include <string>
#include <optional>
#include "MatrixTools.h"
#include <variant>
#include <iostream>

namespace MatrixProject {
/**
 * @brief 计算最小的 n，使得 2^n >= x
 * @param x 输入值
 * @return 最小的 n 值
 */
int smallest_n_for_2n_geq_x(double x);

/**
 * @brief 计算矩阵某一列（不包括对角线）的非零元素个数
 * @param matrix 输入矩阵
 * @param col 列索引
 * @return 非零元素的个数
 */
size_t count_non_zero_elements_in_column(const Matrix& matrix, size_t col);

/**
 * @brief 计算矩阵对角线上某元素的非零元素个数
 * @param matrix 输入矩阵
 * @param i 行索引（也是列索引）
 * @return 非零元素的个数
 */
size_t count_non_zero_elements_in_column_diag(const Matrix& matrix, size_t i);

/**
 * @brief 计算移位次数
 * @param value 输入值
 * @return 移位次数
 */
int calculate_shift(double value);

/**
 * @brief 将矩阵某一列的元素（不包括对角线）添加到列表中
 * @param A 输入矩阵
 * @param i 列索引
 * @return 包含列元素的列表
 */
std::vector<double> add_column_to_list(const Matrix& A, size_t i);

/**
 * @brief 将矩阵某一列的每一行元素添加到列表中
 * @param A 输入矩阵
 * @param i 列索引
 * @param len 列长度
 * @return 包含列元素的列表
 */
std::vector<Element> add_column_to_list_each_row(const Matrix& A, size_t i, size_t len);

/**
 * @brief 将矩阵某一列的对角线元素添加到列表中
 * @param A 输入矩阵
 * @param i 列索引
 * @return 包含对角线元素的列表
 */
std::vector<double> add_column_to_list_diag(const Matrix& A, size_t i);

/**
 * @brief 计算列表中指定索引之前的元素总数
 * @param lst 输入列表
 * @param i 指定索引
 * @return 元素总数
 */
size_t count_elements_before_index(const std::vector<Element>& lst, size_t i);

/**
 * @brief 在每行列表中查找元素的位置
 * @param matrix 输入矩阵
 * @param list_vals 列表值
 * @param list_each_row 每行的列表
 * @param col 列索引
 * @param row 行索引
 * @return 元素的位置（如果找到），否则返回 std::nullopt
 */
std::optional<size_t> find_element_in_list_each_row(
    const Matrix& matrix,
    const std::vector<double>& list_vals,
    const std::vector<Element>& list_each_row,
    size_t col,
    size_t row
);

/**
 * @brief 将整数转换为固定宽度的二进制字符串
 * @param n 整数
 * @param width 位宽
 * @return 二进制字符串
 */
std::string int_to_binary_fixed_width(int n, int width);

/**
 * @brief 对二进制字符串列表进行填充，使其长度一致
 * @param binary_list 输入的二进制字符串列表
 * @return 填充后的二进制字符串列表
 */
std::vector<std::string> pad_binary_list(const std::vector<std::string>& binary_list);

/**
 * @brief 获取矩阵中最后一个非零对角线元素的索引
 * @param matrix 输入矩阵
 * @return 最大非零对角线元素的索引
 */
int max_nonzero_diagonal_index(const Matrix& matrix);

/**
 * @brief 计算矩阵对角线上非零元素的个数
 * @param matrix 输入矩阵
 * @return 非零对角线元素的个数
 */
size_t count_nonzero_diagonal(const Matrix& matrix);

/**
 * @brief 生成选择信号
 * @param matrix_list 矩阵列表
 * @param matrix_a 矩阵 A
 * @param matrix_b 矩阵 B
 * @return 选择信号的字符串列表
 */

// 函数声明
std::vector<Element> add_column_to_list_with_list(const Matrix& A, size_t i);

std::vector<std::variant<int, std::vector<int>>> compute_powers(const std::vector<Element>& shift_list);

int max_input_width(const std::vector<std::variant<int, std::vector<int>>>& new_shift, const std::vector<int>& width_v);

int gcd_list(const std::vector<double>& list);

std::string gen_verilog_expression(const MatrixGroup& matrix_list, const Matrix& matrix_a, size_t i, size_t num_matrix, const std::vector<int>& Co_Shift_list);

std::string gen_verilog_expression_fraction(const MatrixGroup& matrix_list, const Matrix& matrix_a, size_t i, size_t num_matrix, std::vector<int> width_v);

std::string generate_compute_nodes(
    int i,
    const std::vector<int>& width_v,
    const std::vector<int>& Co_Shift_l,
    const std::vector<int>& Co_Shift_r,
    const std::vector<int>& width_l_input_list,
    const std::vector<int>& width_r_input_list
);


std::string generate_compute_nodes_fraction(
    int i, 
    const std::vector<int>& width_v, 
    const std::vector<int>& width_l_input_list,
    const std::vector<int>& width_r_input_list
);

std::vector<int> calculate_y_width(const std::vector<std::string>& output_assignments, const std::vector<int>& v_widths, int num_y);

// Flatten Elements and keep only non-zero values
std::vector<double> flatten_elements(const std::vector<Element>& elements);

// 打印 std::vector<Element> 类型的辅助函数
void print_shift(const std::vector<Element>& shift_l);

// 打印 Element 类型的辅助函数
void print_element(const Element& elem);
} // namespace MatrixProject

#endif // MATRIX_TOOLS_H