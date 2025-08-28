#include "MatrixAB.h"
#include <iostream>
#include <cmath>
#include <algorithm>

namespace MatrixProject {

// 辅助函数：检查 Element 是否为零
bool is_zero(const Element& elem) {
    if (std::holds_alternative<double>(elem)) {
        return std::get<double>(elem) == 0.0;
    } else if (std::holds_alternative<std::vector<double>>(elem)) {
        const auto& vec = std::get<std::vector<double>>(elem);
        return vec.empty() || std::all_of(vec.begin(), vec.end(), [](double val) { return val == 0.0; });
    }
    return true; // 如果没有持有值，视为零
}

// 辅助函数：检查 ElementType 是否等于另一个 ElementType
bool element_equals(const Element& elem1, const Element& elem2) {
    if (std::holds_alternative<double>(elem1) && std::holds_alternative<double>(elem2)) {
        return std::get<double>(elem1) == std::get<double>(elem2);
    } else if (std::holds_alternative<std::vector<double>>(elem1) && std::holds_alternative<std::vector<double>>(elem2)) {
        return std::get<std::vector<double>>(elem1) == std::get<std::vector<double>>(elem2);
    }
    return false;
}

bool element_in_column(const Matrix& matrix, double element, size_t col) {
    for (const auto& row : matrix) {
        const Element& item = row[col];
        if (std::holds_alternative<std::vector<double>>(item)) {
            const auto& vec = std::get<std::vector<double>>(item);
            if (std::find(vec.begin(), vec.end(), element) != vec.end()) {
                return true;
            }
        } else if (std::holds_alternative<double>(item)) {
            if (std::get<double>(item) == element) {
                return true;
            }
        }
    }
    return false;
}

bool element_not_in_position(const Element& position, double element) {
    if (std::holds_alternative<std::vector<double>>(position)) {
        const auto& vec = std::get<std::vector<double>>(position);
        return std::find(vec.begin(), vec.end(), element) == vec.end();
    } else if (std::holds_alternative<double>(position)) {
        return std::get<double>(position) != element;
    }
    return true;
}


std::tuple<Matrix, Matrix> merge_matrices_custom(const std::vector<Matrix>& matrices, bool needPrint) {
    // 找出最大矩阵的尺寸，用于统一尺寸
    size_t max_rows = 0;
    size_t max_cols = 0;
    for (const auto& matrix : matrices) {
        max_rows = std::max(max_rows, matrix.size());
        if (!matrix.empty()) {
            max_cols = std::max(max_cols, matrix[0].size());
        }
    }

    // 初始化两个矩阵 A 和 B，大小为最大尺寸，初始值为 0.0
    Matrix A(max_rows, std::vector<Element>(max_cols, 0.0));
    Matrix B(max_rows, std::vector<Element>(max_cols, 0.0));

    for (const auto& matrix : matrices) {
        size_t rows = matrix.size();
        size_t cols = matrix[0].size();
        //std::cout << "selected matrix :" <<  std::endl;
        //print_matrix(remove_duplicates_from_vectors(matrix));

        for (size_t j = 0; j < cols; ++j) { // 按列处理
            // 找到这一列的两个数字（不一定在同一行）
            std::vector<std::pair<double, size_t>> values;
            for (size_t i = 0; i < rows; ++i) {
                // 如果是对角线上的数字，直接放到矩阵 A 和 B 的相同位置
                if (i == j) {
                    const Element& elem = matrix[i][j];
                    if ((std::holds_alternative<double>(elem) || std::holds_alternative<std::vector<double>>(elem)) && !is_zero(elem)) {
                        // 放入 A[i][j] 和 B[i][j]
                        auto place_in_matrix = [&](Matrix& mat) {
                            if (is_zero(mat[i][j])) {
                                mat[i][j] = elem;
                            } else {
                                if (std::holds_alternative<double>(mat[i][j])) {
                                    std::vector<double> vec = {std::get<double>(mat[i][j])};
                                    if (std::holds_alternative<double>(elem)) {
                                        if (element_not_in_position(mat[i][j], std::get<double>(elem))) {
                                            vec.push_back(std::get<double>(elem));
                                        }
                                    } else {
                                        const auto& elem_vec = std::get<std::vector<double>>(elem);
                                        vec.insert(vec.end(), elem_vec.begin(), elem_vec.end());
                                    }
                                    mat[i][j] = vec;
                                } else {
                                    auto& mat_vec = std::get<std::vector<double>>(mat[i][j]);
                                    if (std::holds_alternative<double>(elem)) {
                                        double val = std::get<double>(elem);
                                        if (element_not_in_position(mat[i][j], val)) {
                                            mat_vec.push_back(val);
                                        }
                                    } else {
                                        const auto& elem_vec = std::get<std::vector<double>>(elem);
                                        mat_vec.insert(mat_vec.end(), elem_vec.begin(), elem_vec.end());
                                    }
                                }
                            }
                        };
                        place_in_matrix(A);
                        place_in_matrix(B);
                    }
                    continue;
                }

                // 检查元素是否是数字或向量，添加到 values
                const Element& elem = matrix[i][j];
                if (std::holds_alternative<double>(elem)) {
                    double val = std::get<double>(elem);
                    if (val != 0.0) {
                        values.emplace_back(val, i);
                    }
                } else if (std::holds_alternative<std::vector<double>>(elem)) {
                    const auto& vec = std::get<std::vector<double>>(elem);
                    if (vec.size() == 2) {
                        values.emplace_back(vec[0], i);
                        values.emplace_back(vec[1], i);
                    } else if (vec.size() == 1 && vec[0] != 0.0) {
                        values.emplace_back(vec[0], i);
                    }
                }
            }

            // 如果这一列没有找到两个数字，跳过
            if (values.size() < 2) {
                continue;
            }

            // 规则 1：如果找到两个负数，报错
            if (values[0].first < 0.0 && values[1].first < 0.0) {
                throw std::runtime_error("第 " + std::to_string(j + 1) + " 列有两个负数，无法处理。");
            }

            // 规则 2 和 3
            double val1 = values[0].first;
            size_t idx1 = values[0].second;
            double val2 = values[1].first;
            size_t idx2 = values[1].second;

            // 检查 A 和 B 中是否包含相同的元素
            bool val1_in_A = element_in_column(A, val1, j);
            bool val2_in_A = element_in_column(A, val2, j);
            bool val1_in_B = element_in_column(B, val1, j);
            bool val2_in_B = element_in_column(B, val2, j);

            auto place_value = [&](Matrix& mat, size_t idx, double val) {
                Element& position = mat[idx][j];
                if (is_zero(position)) {
                    position = val;
                } else {
                    if (std::holds_alternative<double>(position)) {
                        double existing_val = std::get<double>(position);
                        if (existing_val != val) {
                            mat[idx][j] = std::vector<double>{existing_val, val};
                        }
                    } else {
                        auto& vec = std::get<std::vector<double>>(position);
                        if (std::find(vec.begin(), vec.end(), val) == vec.end()) {
                            vec.push_back(val);
                        }
                    }
                }
            };

            if (val1_in_A || val2_in_A) {
                // 如果 A 中有相同的数字
                if (val1_in_A) {
                    place_value(A, idx1, val1);
                    place_value(B, idx2, val2);
                } else {
                    place_value(A, idx2, val2);
                    place_value(B, idx1, val1);
                }
            } else if (val1_in_B || val2_in_B) {
                // 如果 B 中有相同的数字
                if (val1_in_B) {
                    place_value(B, idx1, val1);
                    place_value(A, idx2, val2);
                } else {
                    place_value(B, idx2, val2);
                    place_value(A, idx1, val1);
                }
            } else {
                // 如果两个数字在 A 和 B 中都没有相同的
                // 小的放入 A，大的放入 B
                if (std::abs(val1) < std::abs(val2)) {
                    place_value(A, idx1, val1);
                    place_value(B, idx2, val2);
                } else {
                    place_value(A, idx2, val2);
                    place_value(B, idx1, val1);
                }
            }
        }
    }

    if (needPrint) {
        // 打印结果
        std::cout << "矩阵 A:" << std::endl;
        print_matrix(remove_duplicates_from_vectors(A));
        std::cout << "\n矩阵 B:" << std::endl;
        print_matrix(remove_duplicates_from_vectors(B));
    }

    return std::make_tuple(remove_duplicates_from_vectors(A), remove_duplicates_from_vectors(B));  // 使用 std::make_tuple 返回
}


} // namespace MatrixProject