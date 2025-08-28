#include "MatrixToFraction.h"
#include <algorithm>
#include <cmath>

namespace MatrixProject {

// 工具函数：从Element中提取所有数值，用于寻找最大值
static std::vector<double> extract_all_values(const Element& element) {
    std::vector<double> values;
    if (std::holds_alternative<double>(element)) {
        values.push_back(std::get<double>(element));
    } else {
        const auto& vec = std::get<std::vector<double>>(element);
        values.insert(values.end(), vec.begin(), vec.end());
    }
    return values;
}

// 对Matrix的第col列和第col行进行等价变换（与apply_transform中的操作一致但简化）
// a为缩放因子：
// - 第col列非对角元素都乘以a(相当于列除以1/a)
// - 第col行的所有元素都除以a(相当于行乘以1/a)
static Matrix column_row_transform(const Matrix& matrix, size_t col, double a) {
    size_t rows = matrix.size();
    size_t cols = matrix[0].size();
    Matrix transformed_matrix = matrix; // 深拷贝

    // 第col列(非对角线元素)乘以 a
    for (size_t row = 0; row < rows; ++row) {
        if (row != col) {
            Element& element = transformed_matrix[row][col];
            if (std::holds_alternative<std::vector<double>>(element)) {
                auto& vec = std::get<std::vector<double>>(element);
                for (double &num : vec) {
                    num *= a;
                }
            } else {
                double &val = std::get<double>(element);
                val *= a;
            }
        }
    }

    // 第col行所有元素除以 a
    for (size_t c = 0; c < cols; ++c) {
        Element& element = transformed_matrix[col][c];
        if (std::holds_alternative<std::vector<double>>(element)) {
            auto& vec = std::get<std::vector<double>>(element);
            for (double &num : vec) {
                num /= a;
            }
        } else {
            double &val = std::get<double>(element);
            val /= a;
        }
    }

    return transformed_matrix;
}

Matrix fractionTransform(const Matrix& inputMatrix) {
    if (inputMatrix.empty() || inputMatrix[0].empty()) {
        return inputMatrix; // 空矩阵直接返回
    }

    // 从第二列（index=1）开始处理
    Matrix result = inputMatrix;
    size_t rows = result.size();
    size_t cols = result[0].size();

    for (size_t col = 1; col < cols; ++col) {
        // 找到该列的最大绝对值
        double maxVal = 0.0;
        bool first_flag = true;
        for (size_t row = 0; row < rows; ++row) {
            const Element& element = result[row][col];
            std::vector<double> vals = extract_all_values(element);
            for (double v : vals) {
                double abs_v = std::abs(v);
                if (first_flag) {
                    maxVal = abs_v;
                    first_flag = false;
                } else {
                    if (abs_v > maxVal) {
                        maxVal = abs_v;
                    }
                }
            }
        }

        if (maxVal == 0.0) {
            // 如果该列中所有元素为0或空，无需变换
            continue;
        }

        // a = 1 / maxVal，进行等价变换
        double a = 1.0 / maxVal;
        result = column_row_transform(result, col, a);
    }

    return result;
}

} // namespace MatrixProject