#include "MatrixLoader.h"
#include <fstream>
#include <sstream>
#include <set>
#include <algorithm>

namespace MatrixProject {

// 对矩阵元素进行排序
MatrixGroups sort_matrix_elements(MatrixGroups& matrix_groups) {
    for (auto& group : matrix_groups) {
        for (auto& matrix : group) {
            for (auto& row : matrix) {
                for (auto& element : row) {
                    if (std::holds_alternative<std::vector<double>>(element)) {
                        auto& vec = std::get<std::vector<double>>(element);
                        std::sort(vec.begin(), vec.end());
                    }
                    // 如果元素是 double，则无需排序
                }
            }
        }
    }
    return matrix_groups;
}

MatrixGroups read_matrix_from_txt(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }

    MatrixGroups matrix_groups;
    MatrixGroup current_matrix_group;
    Matrix current_matrix;
    std::vector<Element> current_row;
    std::string line;

    while (std::getline(file, line)) {
        // 去除行首尾空白字符
        line.erase(0, line.find_first_not_of(" \t\n\r"));
        line.erase(line.find_last_not_of(" \t\n\r") + 1);

        if (line == "done") {
            // 结束当前矩阵组
            if (!current_matrix.empty()) {
                current_matrix_group.push_back(current_matrix);
                current_matrix.clear();
            }
            if (!current_matrix_group.empty()) {
                matrix_groups.push_back(current_matrix_group);
                current_matrix_group.clear();
            }
            continue;
        }

        if (line.empty()) {
            // 空行表示当前行结束，将其累加到矩阵
            if (!current_row.empty()) {
                current_matrix.push_back(current_row);
                current_row.clear();

                // 检测是否构成一个方阵
                if (current_matrix.size() == current_matrix[0].size()) {
                    current_matrix_group.push_back(current_matrix);
                    current_matrix.clear();
                }
            }
            continue;
        }

        // 处理数据行
        std::istringstream iss(line);
        std::vector<double> values;
        double value;
        while (iss >> value) {
            values.push_back(value);
        }

        if (!values.empty()) {
            // 将数值列表作为一个元素存储到当前行
            current_row.push_back(Element{values});
        }
    }

    // 处理文件末尾可能残留的数据
    if (!current_row.empty()) {
        current_matrix.push_back(current_row);
        current_row.clear();

        // 如果最后的矩阵是方阵，存储它
        if (current_matrix.size() == current_matrix[0].size()) {
            current_matrix_group.push_back(current_matrix);
            current_matrix.clear();
        }
    }

    if (!current_matrix.empty()) {
        current_matrix_group.push_back(current_matrix);
    }
    if (!current_matrix_group.empty()) {
        matrix_groups.push_back(current_matrix_group);
    }

    // 对矩阵元素进行排序（假设外部定义了这个函数）
    matrix_groups = sort_matrix_elements(matrix_groups);

    return matrix_groups;
}


MatrixGroups remove_duplicate_matrices(const MatrixGroups& matrix_groups) {
    MatrixGroups unique_matrix_groups;
    for (const auto& group : matrix_groups) {
        MatrixGroup unique_matrices;
        std::set<std::string> seen_matrices;
        for (const auto& matrix : group) {
            std::ostringstream oss;
            for (const auto& row : matrix) {
                for (const auto& elem : row) {
                    if (std::holds_alternative<double>(elem)) {
                        oss << std::get<double>(elem) << ",";
                    } else if (std::holds_alternative<std::vector<double>>(elem)) {
                        const auto& vec = std::get<std::vector<double>>(elem);
                        for (const auto& num : vec) {
                            oss << num << ",";
                        }
                    }
                }
                oss << ";";
            }
            std::string matrix_str = oss.str();
            if (seen_matrices.find(matrix_str) == seen_matrices.end()) {
                unique_matrices.push_back(matrix);
                seen_matrices.insert(matrix_str);
            }
        }
        unique_matrix_groups.push_back(unique_matrices);
    }
    return unique_matrix_groups;
}

} // namespace MatrixProject