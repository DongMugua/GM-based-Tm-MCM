#ifndef MATRIX_TRANSFORM_H
#define MATRIX_TRANSFORM_H

#include "MatrixLoader.h"
#include "MatrixTools.h"
#include <vector>

namespace MatrixProject {

    /**
     * @brief 计算 L 值，即 x1 和 x2 中最小绝对值的倒数
     * @param x1 数字 x1
     * @param x2 数字 x2
     * @return 计算得到的 L 值
     */
    double compute_L(double x1, double x2);

    /**
     * @brief 提取矩阵某一列中的两个非对角线非零元素
     * @param matrix 输入矩阵
     * @param col 列索引
     * @param rows 矩阵的行数
     * @return 包含两个非零元素的向量
     */
    std::vector<double> extract_non_diag_elements(const Matrix& matrix, size_t col, size_t rows);

    /**
     * @brief 对单个元素或列表进行除法操作
     * @param element 要处理的元素
     * @param a 除数 a
     * @return 处理后的元素
     */
    Element divide_element_by_a(const Element& element, double a);

    /**
     * @brief 对矩阵执行等价变换
     * @param matrix 输入矩阵
     * @param col 列索引
     * @param a 变换系数 a
     * @param minimum_value 最小值，用于检查元素是否过小
     * @return 一个 pair，第一个元素是变换后的矩阵，第二个元素是布尔值，表示是否达到最小值
     */
    std::pair<Matrix, bool> apply_transform(const Matrix& matrix, size_t col, double a, double minimum_value);

    /**
     * @brief 递归处理矩阵的列，进行变换
     * @param matrix 输入矩阵
     * @param current_col 当前处理的列索引
     * @param cols 总列数
     * @param matrix_group 存储生成的矩阵组
     * @param minimum_value 最小值，用于检查元素是否过小
     */
    void process_columns(const Matrix& matrix, size_t current_col, size_t cols, std::vector<Matrix>& matrix_group, double minimum_value);

    /**
     * @brief 对矩阵组执行变换操作
     * @param matrix_group 输入的矩阵组
     * @param minimum_value 最小值，用于检查元素是否过小
     * @return 变换后的矩阵组
     */
    std::vector<Matrix> matrix_transform(const std::vector<Matrix>& matrix_group, double minimum_value);

} // namespace MatrixProject

#endif // MATRIX_TRANSFORM_H