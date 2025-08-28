#include "MatrixTransform.h"
#include <algorithm>
#include <cmath>
#include <iostream>

namespace MatrixProject {

double compute_L(double x1, double x2) {
    if (x1 == 0.0 || x2 == 0.0) {
        return 1.0; // 避免除以 0
    }

    double min_abs = std::min(std::abs(x1), std::abs(x2));
    if (min_abs >= 1.0) {
        return 1.0 / min_abs; // L 为 x1 和 x2 中最小值的倒数
    } else {
        return 1.0;
    }
}

std::vector<double> extract_non_diag_elements(const Matrix& matrix, size_t col, size_t rows) {
    std::vector<double> non_diag_elements;
    for (size_t row = 0; row < rows; ++row) {
        if (row == col) {
            continue;
        }
        const Element& element = matrix[row][col];

        // 检查是否为向量，并且元素不为零
        if (std::holds_alternative<std::vector<double>>(element)) {
            const auto& vec = std::get<std::vector<double>>(element);
            if (vec.size() == 2) {
                for (double e : vec) {
                    if (e != 0.0) {
                        non_diag_elements.push_back(e);
                    }
                }
            } else if (vec.size() == 1 && vec[0] != 0.0) {
                non_diag_elements.push_back(vec[0]);
            }
        } else if (std::holds_alternative<double>(element)) {
            double val = std::get<double>(element);
            if (val != 0.0) {
                non_diag_elements.push_back(val);
            }
        }
    }

    // 只保留两个非零数字
    if (non_diag_elements.size() > 2) {
        non_diag_elements.resize(2);
    }

    return non_diag_elements;
}

Element divide_element_by_a(const Element& element, double a) {
    if (std::holds_alternative<std::vector<double>>(element)) {
        const auto& vec = std::get<std::vector<double>>(element);
        std::vector<double> result;
        for (double num : vec) {
            result.push_back(num / a);
        }
        return result;
    } else if (std::holds_alternative<double>(element)) {
        double val = std::get<double>(element);
        return val / a;
    }
    return element;
}

std::pair<Matrix, bool> apply_transform(const Matrix& matrix, size_t col, double a, double minimum_value) {
    size_t rows = matrix.size();
    Matrix transformed_matrix = matrix; // 深拷贝
    bool reach_minimum = false;

    // 对第 col 列的非对角线元素乘以 a
    for (size_t row = 0; row < rows; ++row) {
        if (row != col) {
            Element& element = transformed_matrix[row][col];

            if (std::holds_alternative<std::vector<double>>(element)) {
                auto& vec = std::get<std::vector<double>>(element);
                for (double& num : vec) {
                    num *= a;
                }
            } else if (std::holds_alternative<double>(element)) {
                double& val = std::get<double>(element);
                val *= a;
            }
        }
    }

    // 对第 col 行的所有元素除以 a
    for (size_t i = 0; i < transformed_matrix[col].size(); ++i) {
        transformed_matrix[col][i] = divide_element_by_a(transformed_matrix[col][i], a);
    }

    // 遍历整个矩阵，检查是否有小于 minimum_value 的元素
    for (size_t row = 0; row < rows; ++row) {
        for (size_t c = 0; c < transformed_matrix[row].size(); ++c) {
            const Element& element = transformed_matrix[row][c];

            if (std::holds_alternative<double>(element)) {
                double val = std::get<double>(element);
                if (std::abs(val) < minimum_value && val != 0.0) {
                    reach_minimum = true;
                    return {transformed_matrix, reach_minimum};
                }
            } else if (std::holds_alternative<std::vector<double>>(element)) {
                const auto& vec = std::get<std::vector<double>>(element);
                for (double num : vec) {
                    if (std::abs(num) < minimum_value && num != 0.0) {
                        reach_minimum = true;
                        return {transformed_matrix, reach_minimum};
                    }
                }
            }
        }
    }

    return {transformed_matrix, reach_minimum};
}

void process_columns(const Matrix& matrix, size_t current_col, size_t cols, std::vector<Matrix>& matrix_group, double minimum_value) {
    if (current_col >= cols) { // 递归终止条件
        return;
    }

    size_t rows = matrix.size();

    // 提取当前列的两个非对角线元素
    std::vector<double> non_diag_elements = extract_non_diag_elements(matrix, current_col, rows);

    if (non_diag_elements.size() != 2) {
        // 如果没有足够的非对角线元素，继续处理下一列
        process_columns(matrix, current_col + 1, cols, matrix_group, minimum_value);
        return;
    }

    double x1 = non_diag_elements[0];
    double x2 = non_diag_elements[1];
    double L = compute_L(x1, x2); // L 为 x1 和 x2 中最小值的倒数

    // 生成 a_values 列表
    std::vector<double> a_values;
    for (int i = 0; i <= 10; ++i) {
        double val = L * std::pow(2, i);
        if (val <= 64.0) {
            a_values.push_back(val);
        } else {
            break;
        }
    }

    // 对当前列的每个 a 值进行处理，并递归处理后续列
    for (double a : a_values) {
        if (current_col < cols) {
            auto [transformed_matrix, reach_minimum] = apply_transform(matrix, current_col, a, minimum_value);
            if (!reach_minimum) {
                matrix_group.push_back(transformed_matrix);
                // 递归处理下一列
                process_columns(transformed_matrix, current_col + 1, cols, matrix_group, minimum_value);
            }
            // if (reach_minimum) {
            //     break;
            // }
        }
    }
}

std::vector<Matrix> matrix_transform(const std::vector<Matrix>& matrix_group, double minimum_value) {
    size_t count_total = 0;
    std::vector<Matrix> matrix_group_result; // 用于保存所有新生成的矩阵

    // 遍历原始矩阵组
    for (const Matrix& matrix : matrix_group) {
        size_t rows = matrix.size();
        size_t cols = matrix[0].size();
        std::vector<Matrix> result_for_current_matrix; // 临时存储每个原始矩阵的生成结果

        // 从第 1 列开始递归处理
        process_columns(matrix, 1, cols, result_for_current_matrix, minimum_value);

        // 将当前矩阵的结果添加到最终结果列表中
        matrix_group_result.insert(matrix_group_result.end(), result_for_current_matrix.begin(), result_for_current_matrix.end());
        count_total += result_for_current_matrix.size();
    }

    std::cout << "新增矩阵个数：" << count_total << std::endl;
    return matrix_group_result;
}

} // namespace MatrixProject