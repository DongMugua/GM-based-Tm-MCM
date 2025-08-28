#include "MatrixTools.h"
#include <iostream>
#include <vector>
#include <variant>
#include <iomanip>
#include <sstream>
#include <set>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <stdexcept>

namespace MatrixProject {

    // 计算单个元素的字符串长度
    size_t element_width(const Element& elem) {
        if (std::holds_alternative<double>(elem)) {
            std::ostringstream oss;
            oss << std::get<double>(elem);
            return oss.str().length();
        } else if (std::holds_alternative<std::vector<double>>(elem)) {
            const auto& vec = std::get<std::vector<double>>(elem);
            std::ostringstream oss;
            oss << "[";
            for (size_t i = 0; i < vec.size(); ++i) {
                oss << vec[i];
                if (i < vec.size() - 1) oss << " ";
            }
            oss << "]";
            return oss.str().length();
        }
        return 0;
    }

    void print_matrix(const MatrixGroup& group) {
        if (group.empty()) return;

        for (std::size_t idx = 0; idx < group.size(); ++idx) {
            std::cout << "Matrix " << idx << ":\n";
            print_matrix(group[idx]);          // reuse the single‑matrix printer
        }
    }


    void print_matrix(const MatrixGroups& groups) {
        if (groups.empty()) return;

        for (std::size_t g = 0; g < groups.size(); ++g) {
            std::cout << "=== Group " << g << " ===\n";
            print_matrix(groups[g]);   // reuse MatrixGroup printer
        }

        std::cout << "====== END ======\n";

    }

    // 打印矩阵
    void print_matrix(const Matrix& matrix) {
        if (matrix.empty()) return;

        // 计算每列最大宽度
        std::vector<size_t> col_widths(matrix[0].size(), 0);
        for (const auto& row : matrix) {
            for (size_t j = 0; j < row.size(); ++j) {
                col_widths[j] = std::max(col_widths[j], element_width(row[j]));
            }
        }

        // 打印矩阵
        for (const auto& row : matrix) {
            for (size_t j = 0; j < row.size(); ++j) {
                const auto& elem = row[j];
                if (std::holds_alternative<double>(elem)) {
                    std::cout << std::setw(col_widths[j]) << std::get<double>(elem) << " ";
                } else if (std::holds_alternative<std::vector<double>>(elem)) {
                    const auto& vec = std::get<std::vector<double>>(elem);
                    std::ostringstream oss;
                    oss << "[";
                    for (size_t i = 0; i < vec.size(); ++i) {
                        oss << vec[i];
                        if (i < vec.size() - 1) oss << " ";
                    }
                    oss << "]";
                    std::cout << std::setw(col_widths[j]) << oss.str() << " ";
                }
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }


// 打印矩阵
void print_matrix(const Matrix& matrix, std::ostream& out) {
    if (matrix.empty()) return;

    // 计算每列最大宽度
    std::vector<size_t> col_widths(matrix[0].size(), 0);
    for (const auto& row : matrix) {
        for (size_t j = 0; j < row.size(); ++j) {
            col_widths[j] = std::max(col_widths[j], element_width(row[j]));
        }
    }

    // 打印矩阵
    for (const auto& row : matrix) {
        for (size_t j = 0; j < row.size(); ++j) {
            const auto& elem = row[j];
            if (std::holds_alternative<double>(elem)) {
                out << std::setw(col_widths[j]) << std::get<double>(elem) << " ";
            } else if (std::holds_alternative<std::vector<double>>(elem)) {
                const auto& vec = std::get<std::vector<double>>(elem);
                std::ostringstream oss;
                oss << "[";
                for (size_t i = 0; i < vec.size(); ++i) {
                    oss << vec[i];
                    if (i < vec.size() - 1) oss << " ";
                }
                oss << "]";
                out << std::setw(col_widths[j]) << oss.str() << " ";
            }
        }
        out << std::endl;
    }
}

    Matrix remove_duplicates_from_vectors(const Matrix& matrix) {
        Matrix new_matrix = matrix; // 复制输入矩阵
        for (auto& row : new_matrix) {
            for (auto& elem : row) {
                if (std::holds_alternative<std::vector<double>>(elem)) {
                    auto& vec = std::get<std::vector<double>>(elem); // 获取向量引用

                    // 判断向量是否包含非零元素
                    bool has_non_zero = std::any_of(vec.begin(), vec.end(), [](double v) { return v != 0.0; });

                    if (has_non_zero) {
                        // 存在非零元素：移除零，并去重
                        std::set<double> unique_elements;
                        for (double val : vec) {
                            if (val != 0.0) {
                                unique_elements.insert(val);
                            }
                        }
                        vec.assign(unique_elements.begin(), unique_elements.end());
                    } else {
                        // 全为零：保留一个零
                        vec = {0.0};
                    }
                }
            }
        }
        return new_matrix;
    }

Matrix ensure_list_format(const Matrix& matrix) {
    Matrix new_matrix = matrix;
    for (auto& row : new_matrix) {
        for (auto& elem : row) {
            if (std::holds_alternative<double>(elem)) {
                double num = std::get<double>(elem);
                elem = std::vector<double>{num};
            }
        }
    }
    return new_matrix;
}

MatrixGroups pad_all_groups(const MatrixGroups& matrix_groups) {
    int max_size = get_max_square_matrix_size(matrix_groups);
    MatrixGroups padded_groups;

    for (const auto& group : matrix_groups) {
        MatrixGroup padded_group;

        for (const auto& matrix : group) {
            if (matrix.size() < max_size) {
                // 使用 iterative_pad_matrix 补全矩阵
                std::vector<Matrix> padded_versions = iterative_pad_matrix(matrix, max_size);
                padded_group.insert(padded_group.end(), padded_versions.begin(), padded_versions.end());
            } else {
                padded_group.push_back(matrix); // 保持最大阶矩阵不变
            }
        }

        padded_groups.push_back(padded_group);
    }

    return padded_groups;
}

int get_max_square_matrix_size(const MatrixGroups& matrix_groups) {
    int max_size = 0;
    for (const auto& group : matrix_groups) {
        for (const auto& matrix : group) {
            if (matrix.size() == matrix[0].size()) {
                max_size = std::max(max_size, static_cast<int>(matrix.size()));
            }
        }
    }
    return max_size;
}

Matrix pad_square_matrix(const Matrix& matrix, int size) {
    Matrix padded_matrix = matrix;
    int current_size = static_cast<int>(matrix.size());
    if (current_size < size) {
        // 添加行
        for (int i = current_size; i < size; ++i) {
            std::vector<Element> new_row(size, 0.0);
            padded_matrix.push_back(new_row);
        }
        // 添加列
        for (auto& row : padded_matrix) {
            for (int j = static_cast<int>(row.size()); j < size; ++j) {
                row.push_back(0.0);
            }
        }
    }
    return padded_matrix;
}

Matrix merge_matrices(const Matrix& matrixA, const Matrix& matrixB) {
    // 检查矩阵是否为空
    if (matrixA.empty() || matrixB.empty()) {
        throw std::runtime_error("One or both matrices are empty");
    }

    // 检查矩阵行数是否相同
    if (matrixA.size() != matrixB.size()) {
        throw std::runtime_error("Matrices must have the same number of rows");
    }

    // 检查每一行的列数是否相同
    for (size_t i = 0; i < matrixA.size(); ++i) {
        if (matrixA[i].size() != matrixB[i].size()) {
            throw std::runtime_error("Matrices must have the same number of columns in each row");
        }
    }

    size_t rows = matrixA.size();
    size_t cols = matrixA[0].size();
    Matrix merged_matrix(rows, std::vector<Element>(cols));

    for (size_t i = 0; i < rows; ++i) {
        cols = matrixA[i].size(); // 更新列数
        for (size_t j = 0; j < cols; ++j) {
            const Element& elemA = matrixA[i][j];
            const Element& elemB = matrixB[i][j];

            std::vector<double> merged_values;

            // 处理 elemA
            if (std::holds_alternative<double>(elemA)) {
                merged_values.push_back(std::get<double>(elemA));
            } else if (std::holds_alternative<std::vector<double>>(elemA)) {
                const auto& vecA = std::get<std::vector<double>>(elemA);
                merged_values.insert(merged_values.end(), vecA.begin(), vecA.end());
            }

            // 处理 elemB
            if (std::holds_alternative<double>(elemB)) {
                merged_values.push_back(std::get<double>(elemB));
            } else if (std::holds_alternative<std::vector<double>>(elemB)) {
                const auto& vecB = std::get<std::vector<double>>(elemB);
                merged_values.insert(merged_values.end(), vecB.begin(), vecB.end());
            }

            // 存入合并矩阵
            merged_matrix[i][j] = merged_values;
        }
    }

    return remove_duplicates_from_vectors(merged_matrix);
}

    std::vector<Matrix> iterative_pad_matrix(const Matrix& matrix, int max_size) {
        std::vector<Matrix> padded_matrices; // 最终结果矩阵列表
        padded_matrices.push_back(matrix);  // 初始矩阵作为第一个版本

        std::vector<Matrix> current_matrices{matrix}; // 当前需要处理的矩阵列表

        // 逐步扩展矩阵，直到达到目标大小
        while (current_matrices[0].size() < max_size) {
            std::vector<Matrix> next_round_matrices;

            for (const auto& mat : current_matrices) {
                int current_size = static_cast<int>(mat.size());

                // 遍历每一列和每一行
                for (int i = 0; i <= current_size; ++i) { // 注意：从索引 0 开始插入
                    Matrix padded_matrix = mat;

                    // 在列 i 插入一列零
                    for (auto& row : padded_matrix) {
                        row.insert(row.begin() + i, 0.0); // 插入零列
                    }

                    // 在行 i 插入一行零
                    std::vector<Element> new_row(padded_matrix[0].size(), 0.0);
                    padded_matrix.insert(padded_matrix.begin() + i, new_row); // 插入零行

                    // 将补全的矩阵加入下一轮处理列表和最终结果
                    next_round_matrices.push_back(padded_matrix);
                    padded_matrices.push_back(padded_matrix);
                }
            }

            // 更新当前矩阵列表
            current_matrices = next_round_matrices;
        }

        return padded_matrices;
    }

    // 检查两个 double 是否相等（考虑容差）
    bool is_close(double a, double b, double tol = 1e-8) {
        return std::fabs(a - b) <= tol;
    }

    // 检查两个 Element 是否相等
    bool elements_equal(const Element& e1, const Element& e2, double tol = 1e-8) {
        if (std::holds_alternative<double>(e1) && std::holds_alternative<double>(e2)) {
            return is_close(std::get<double>(e1), std::get<double>(e2), tol); // 浮点数比较
        }
        if (std::holds_alternative<std::vector<double>>(e1) && std::holds_alternative<std::vector<double>>(e2)) {
            auto vec1 = std::get<std::vector<double>>(e1);
            auto vec2 = std::get<std::vector<double>>(e2);
            std::sort(vec1.begin(), vec1.end()); // 对向量排序
            std::sort(vec2.begin(), vec2.end());
            return vec1 == vec2; // 排序后比较
        }
        return false; // 类型不匹配
    }

    // 检查矩阵是否存在于矩阵组中
    bool matrix_exists(const Matrix& matrix, const MatrixGroup& matrix_group, double tol = 1e-8) {
        for (const auto& mat : matrix_group) {
            if (matrix.size() != mat.size() || matrix[0].size() != mat[0].size()) {
                continue; // 跳过大小不匹配的矩阵
            }

            bool match = true;
            for (size_t i = 0; i < matrix.size(); ++i) {
                for (size_t j = 0; j < matrix[i].size(); ++j) {
                    if (!elements_equal(matrix[i][j], mat[i][j], tol)) {
                        match = false; // 如果元素不匹配，标记为不相等
                        break;
                    }
                }
                if (!match) break;
            }

            if (match) return true; // 如果找到匹配的矩阵，直接返回
        }

        return false; // 未找到匹配的矩阵
    }

} // namespace MatrixProject