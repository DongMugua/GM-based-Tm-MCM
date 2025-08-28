#include "MatrixCost.h"
#include "MatrixTools.h"
#include "Config.h"

#include <cmath>
#include <set>
#include <map>
#include <numeric>
#include <iostream>
#include <algorithm>

namespace MatrixProject {

std::vector<int> count_non_zero_elements_off_diagonal(const Matrix& matrix) {
    size_t rows = matrix.size();
    size_t cols = matrix[0].size();
    std::vector<int> non_zero_count_vector;

    for (size_t col = 0; col < cols; ++col) {
        std::map<double, int> unique_non_zero_numbers;

        for (size_t row = 0; row < rows; ++row) {
            if (row != col) {  // 排除对角线元素
                const Element& elem = matrix[row][col];

                if (std::holds_alternative<double>(elem)) {
                    double value = std::get<double>(elem);
                    if (value != 0.0) {
                        unique_non_zero_numbers[std::abs(value)] += 1;
                    }
                } else if (std::holds_alternative<std::vector<double>>(elem)) {
                    const auto& vec = std::get<std::vector<double>>(elem);
                    for (double num : vec) {
                        if (num != 0.0) {
                            unique_non_zero_numbers[std::abs(num)] += 1;
                        }
                    }
                }
            }
        }

        // 统计非零数字的总出现次数
        int non_zero_count = std::accumulate(unique_non_zero_numbers.begin(), unique_non_zero_numbers.end(), 0,
                                             [](int sum, const std::pair<double, int>& p) { return sum + p.second; });

        // 如果非零数字个数为1，则将其置为0
        non_zero_count_vector.push_back(non_zero_count == 1 ? 0 : non_zero_count);
    }

    return non_zero_count_vector;
}

double get_min_value(const Element& element, bool return_abs) {
    if (std::holds_alternative<double>(element)) {
        // 如果是单个 double
        double value = std::get<double>(element);
        return return_abs ? std::abs(value) : value;
    } else if (std::holds_alternative<std::vector<double>>(element)) {
        // 如果是一个 vector<double>
        const auto& vec = std::get<std::vector<double>>(element);
        if (vec.empty()) return 0.0;

        // 使用 lambda 比较器基于绝对值进行最小比较
        auto min_it = std::min_element(
            vec.begin(), vec.end(),
            [](double a, double b) {
                return std::abs(a) < std::abs(b);
            }
        );
        // 根据 return_abs 决定返回绝对值还是原值
        return return_abs ? std::abs(*min_it) : *min_it;
    } else {
        throw std::runtime_error("Unsupported element type in get_min_value");
    }
}

double get_max_value(const Element& element, bool return_abs) {
    if (std::holds_alternative<double>(element)) {
        double value = std::get<double>(element);
        return return_abs ? std::abs(value) : value;
    } else if (std::holds_alternative<std::vector<double>>(element)) {
        const auto& vec = std::get<std::vector<double>>(element);
        if (vec.empty()) return 0.0;
        auto max_it = std::max_element(vec.begin(), vec.end(), [](double a, double b) { return std::abs(a) < std::abs(b); });
        return return_abs ? std::abs(*max_it) : *max_it;
    } else {
        throw std::runtime_error("Unsupported element type in get_max_value");
    }
}

Matrix transform_matrix(const Matrix& matrix) {
    Matrix transformed_matrix = matrix;  // 深拷贝
    size_t rows = transformed_matrix.size();
    size_t cols = transformed_matrix[0].size();
    for (size_t col = 1; col < cols - 1; ++col) {
        double fractor_value;
        if (!config.FRACTION_INPUT) {
            fractor_value = 0.0;
            // 获取当前列中所有对角线以上部分的绝对值最大值
            for (size_t row = 0; row < col; ++row) {
                double val = get_max_value(transformed_matrix[row][col], true);
                if (val > fractor_value) {
                    fractor_value = val;
                }
            }
        } else {
            fractor_value = 1;
            // 获取当前列中所有对角线以上部分的绝对值最小值
            for (size_t row = 0; row < col; ++row) {
                double val = get_min_value(transformed_matrix[row][col], true);
                if (val < fractor_value && val != 0.0) {
                    fractor_value = val;
                }
            }
        }
        // 将当前列的绝对值最大值与后面每一列的第 i 行除对角线外的所有值相乘
        for (size_t next_col = col + 1; next_col < cols; ++next_col) {
            size_t i = col;
            if (i == next_col) continue;

            Element& elem = transformed_matrix[i][next_col];
            if (std::holds_alternative<double>(elem)) {
                double value = std::get<double>(elem);
                value *= fractor_value;
                elem = value;
            } else if (std::holds_alternative<std::vector<double>>(elem)) {
                auto& vec = std::get<std::vector<double>>(elem);
                for (double& e : vec) {
                    e *= fractor_value;
                }
            }
        }
    }

    return transformed_matrix;
}

std::vector<double> get_column_vector(const Matrix& matrix) {
    size_t rows = matrix.size();
    size_t cols = matrix[0].size();
    std::vector<double> column_vector(cols, -std::numeric_limits<double>::infinity());

    if (!config.FRACTION_INPUT) {
        for (size_t col = 0; col < cols; ++col) {
            double col_max = -std::numeric_limits<double>::infinity();
            for (size_t row = 0; row < rows; ++row) {
                if (row == col) continue;
                double element_max = get_max_value(matrix[row][col]);
                if (element_max > col_max) {
                    col_max = element_max;
                }
            }

            // if (col_max == -std::numeric_limits<double>::infinity())
            //     col_max = 1;

            column_vector[col] = col_max;
        }
    } else {
        for (size_t col = 0; col < cols; ++col) {
            double col_min = std::numeric_limits<double>::infinity();
            for (size_t row = 0; row < rows; ++row) {
                if (row == col) continue;
                double element_min = get_min_value(matrix[row][col],1);
                if (element_min < col_min && element_min != 0.0) {
                    col_min = element_min;
                }
            }
            if (col_min == std::numeric_limits<double>::infinity())
                col_min = 1;

            column_vector[col] = col_min;
        }
    }

    return column_vector;
}

std::pair<double, double> cost_matrices_mux(const Matrix& matrixA, const Matrix& matrixB, const std::vector<double>& mac_vector) {
    // 由于矩阵可能是稀疏的，需要小心处理元素
    size_t rows = matrixA.size();
    size_t cols = matrixA[0].size();

    // 初始化结果矩阵
    Matrix result_matrix_A(rows, std::vector<Element>(cols));
    Matrix result_matrix_B(rows, std::vector<Element>(cols));

    // 定义一个辅助函数，用于处理每个元素的乘法
    auto multiply_element = [](const Element& element, double multiplier) -> Element {
        if (std::holds_alternative<double>(element)) {
            double value = std::get<double>(element);
            return value * multiplier;
        } else if (std::holds_alternative<std::vector<double>>(element)) {
            const auto& vec = std::get<std::vector<double>>(element);
            std::vector<double> result;
            for (double e : vec) {
                result.push_back(e * multiplier);
            }
            return result;
        }
        return 0.0;
    };

    // 逐列相乘
    for (size_t col = 0; col < cols; ++col) {
        for (size_t row = 0; row < rows; ++row) {
            result_matrix_A[row][col] = multiply_element(matrixA[row][col], mac_vector[col]);
            result_matrix_B[row][col] = multiply_element(matrixB[row][col], mac_vector[col]);
        }
    }

    auto vectorA = get_column_vector(result_matrix_A);
    auto vectorB = get_column_vector(result_matrix_B);

    auto resultA = count_non_zero_elements_off_diagonal(result_matrix_A);
    auto resultB = count_non_zero_elements_off_diagonal(result_matrix_B);

    std::vector<double> MUX_cost_vectorA;
    std::vector<double> MUX_cost_vectorB;

    for (size_t i = 0; i < vectorA.size(); ++i) {
        double mac = vectorA[i];
        int cost = resultA[i];
        double mux_cost = (mac != 0.0) ? cost * (abs(std::log2(mac)) + config.BIT_WIDTH) : 0.0;
        MUX_cost_vectorA.push_back(mux_cost);
    }

    for (size_t i = 0; i < vectorB.size(); ++i) {
        double mac = vectorB[i];
        int cost = resultB[i];
        double mux_cost = (mac != 0.0) ? cost * (abs(std::log2(mac)) + config.BIT_WIDTH) : 0.0;
        MUX_cost_vectorB.push_back(mux_cost);
    }

    // 计算对角线上的非零元素
    int non_zero_count = 0;
    double max_value = 0.0;
    for (size_t i = 0; i < std::min(rows, cols); ++i) {
        const Element& elem = result_matrix_A[i][i];
        if (std::holds_alternative<double>(elem)) {
            double value = std::get<double>(elem);
            if (value != 0.0) {
                non_zero_count += 1;
                if (value > max_value) {
                    max_value = value;
                }
            }
        } else if (std::holds_alternative<std::vector<double>>(elem)) {
            const auto& vec = std::get<std::vector<double>>(elem);
            for (double num : vec) {
                if (num != 0.0) {
                    non_zero_count += 1;
                    if (num > max_value) {
                        max_value = num;
                    }
                }
            }
        }
    }

    if (non_zero_count == 1) {
        non_zero_count = 0;
    }

    // 计算 MUX 数量
    auto count_mux_num = [](const std::vector<int>& array1, const std::vector<int>& array2) -> double {
        double count = 0.0;
        for (int x : array1) {
            if (x != 0 && x != 1) {
                count += x / 2.0;
            }
        }
        for (int x : array2) {
            if (x != 0 && x != 1) {
                count += x / 2.0;
            }
        }
        return count;
    };

    double MUX_num = count_mux_num(resultA, resultB) + ((non_zero_count != 0) ? non_zero_count / 2.0 : 0.0);

    double total_cost = (std::accumulate(MUX_cost_vectorA.begin(), MUX_cost_vectorA.end(), 0.0) +
                         std::accumulate(MUX_cost_vectorB.begin(), MUX_cost_vectorB.end(), 0.0)) * 14.0;

    // TODO 简单地将对角线上的元素按照最大的值来计算的，先注释掉，后面再详细优化
    // if (max_value != 0.0) {
    //     total_cost += non_zero_count * (std::log2(max_value) + config.BIT_WIDTH) * 14.0;
    // }

    return std::make_pair(total_cost, MUX_num);
}

std::vector<int> cost_matrices_adders(const Matrix& matrixA, const Matrix& matrixB) {
    // 确保矩阵大小相同
    if (matrixA.size() != matrixB.size() || matrixA[0].size() != matrixB[0].size()) {
        throw std::runtime_error("Matrices must be of the same size");
    }

    size_t n = matrixA.size();
    std::vector<int> result_vector;

    // 定义辅助函数提取元素中的数字
    auto extract_numbers = [](const Element& element) -> std::vector<double> {
        std::vector<double> nums;
        if (std::holds_alternative<double>(element)) {
            nums.push_back(std::get<double>(element));
        } else if (std::holds_alternative<std::vector<double>>(element)) {
            nums = std::get<std::vector<double>>(element);
        }
        return nums;
    };

    for (size_t col = 0; col < n; ++col) {
        std::vector<double> a_numbers;
        std::vector<double> b_numbers;

        for (size_t row = 0; row < col; ++row) {
            auto nums_a = extract_numbers(matrixA[row][col]);
            auto nums_b = extract_numbers(matrixB[row][col]);
            a_numbers.insert(a_numbers.end(), nums_a.begin(), nums_a.end());
            b_numbers.insert(b_numbers.end(), nums_b.begin(), nums_b.end());
        }

        // 移除零值
        a_numbers.erase(std::remove(a_numbers.begin(), a_numbers.end(), 0.0), a_numbers.end());
        b_numbers.erase(std::remove(b_numbers.begin(), b_numbers.end(), 0.0), b_numbers.end());

        // 判定条件
        bool all_zero = a_numbers.empty() && b_numbers.empty();
        bool all_positive = std::all_of(a_numbers.begin(), a_numbers.end(), [](double num) { return num > 0.0; }) &&
                            std::all_of(b_numbers.begin(), b_numbers.end(), [](double num) { return num > 0.0; });
        bool mix_sign = (std::all_of(a_numbers.begin(), a_numbers.end(), [](double num) { return num > 0.0; }) &&
                         std::all_of(b_numbers.begin(), b_numbers.end(), [](double num) { return num < 0.0; })) ||
                        (std::all_of(a_numbers.begin(), a_numbers.end(), [](double num) { return num < 0.0; }) &&
                         std::all_of(b_numbers.begin(), b_numbers.end(), [](double num) { return num > 0.0; }));

        if (all_zero) {
            result_vector.push_back(0);
        } else if (all_positive) {
            result_vector.push_back(67);  // ADDER COST
        } else if (mix_sign) {
            result_vector.push_back(75);  // SUBTRACTOR COST
        } else {
            result_vector.push_back(93);  // ADD/SUB COST
        }
    }

    return result_vector;
}

std::pair<double, double> cost_analysis(const Matrix& matrixA, const Matrix& matrixB) {
    Matrix union_matrix = merge_matrices(matrixA, matrixB);
    Matrix transformed_matrix = transform_matrix(union_matrix);
    auto mac_vector = get_column_vector(transformed_matrix);
    auto adder_concept_vector = cost_matrices_adders(matrixA, matrixB);

    std::vector<double> adder_cost_vector;
    for (size_t i = 0; i < adder_concept_vector.size(); ++i) {
        int cost = adder_concept_vector[i];
        double mac = mac_vector[i];
        double adder_cost = 0.0;
        if (mac != 0.0) {
            // TODO 后续可以细化第一级加法/减法器的位宽与移位的关系
            if (i == 1 && cost == 67) { // 只有在第一个加法器且为加法时，可以直接近似成bitWidth
                adder_cost = cost * (std::log2(1.0) + config.BIT_WIDTH);
            } else {
                adder_cost = cost * (abs(std::log2(mac)) + config.BIT_WIDTH);
            }
        }

        adder_cost_vector.push_back(adder_cost);
    }

    double total_adder_cost = std::accumulate(adder_cost_vector.begin(), adder_cost_vector.end(), 0.0);
    auto mux_cost = cost_matrices_mux(matrixA, matrixB, mac_vector);
    // double total_cost = (total_adder_cost + mux_cost.first) * mux_cost.second;

    double total_cost = total_adder_cost + mux_cost.first;

    return std::make_pair(total_cost, mux_cost.second);
}

} // namespace MatrixProject