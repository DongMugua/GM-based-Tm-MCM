#include "MatrixToHDL.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cctype>
#include <stdexcept>

namespace MatrixProject {

    // ========== 辅助函数 ==========

    // 去除字符串首尾空白字符
    static std::string trim(const std::string& s) {
        size_t start = s.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) return "";
        size_t end = s.find_last_not_of(" \t\r\n");
        return s.substr(start, end - start + 1);
    }

    // 按空白字符分割字符串
    static std::vector<std::string> split(const std::string& s) {
        std::istringstream iss(s);
        std::vector<std::string> tokens;
        std::string token;
        while (iss >> token) {
            tokens.push_back(token);
        }
        return tokens;
    }

    // 将一行字符串解析为多个 Element
    static std::vector<Element> parseRow(const std::string& line) {
        std::vector<Element> row;
        std::vector<std::string> tokens;
        bool inBracket = false;
        std::string token;

        // 按字符扫描，区分方括号内外
        for (char ch : line) {
            if (ch == '[') {
                // 若不在方括号内且 token 不为空，则先保存
                if (!token.empty() && !inBracket) {
                    tokens.push_back(token);
                    token.clear();
                }
                inBracket = true;
                token.push_back(ch);
            } else if (ch == ']') {
                token.push_back(ch);
                tokens.push_back(token);
                token.clear();
                inBracket = false;
            } else if (std::isspace(static_cast<unsigned char>(ch)) && !inBracket) {
                if (!token.empty()) {
                    tokens.push_back(token);
                    token.clear();
                }
            } else {
                token.push_back(ch);
            }
        }
        if (!token.empty()) {
            tokens.push_back(token);
        }

        // 解析 tokens
        for (auto& tk : tokens) {
            tk = trim(tk);
            if (tk.empty()) continue;
            // 如果以 '[' 开头且以 ']' 结尾，则解析为 vector<double>
            if (tk.front() == '[' && tk.back() == ']') {
                std::string inner = tk.substr(1, tk.size() - 2); // 去掉方括号
                std::vector<std::string> numStrs = split(inner);
                std::vector<double> vec;
                for (const auto& numStr : numStrs) {
                    try {
                        vec.push_back(std::stod(numStr));
                    } catch (...) {
                        std::cerr << "解析数字失败: " << numStr << std::endl;
                    }
                }
                row.push_back(vec);
            } else {
                // 否则当作单个 double
                try {
                    double val = std::stod(tk);
                    row.push_back(val);
                } catch (...) {
                    std::cerr << "解析数字失败: " << tk << std::endl;
                }
            }
        }
        return row;
    }

    bool parse_matrices_from_file(const std::string& filename,
                                  std::vector<Matrix>& selected_matrices,
                                  Matrix& matrixA,
                                  Matrix& matrixB)
    {
        std::ifstream infile(filename);
        if (!infile.is_open()) {
            std::cerr << "无法打开文件: " << filename << std::endl;
            return false;
        }

        // 清空输出参数
        selected_matrices.clear();
        matrixA.clear();
        matrixB.clear();

        // 用于暂存 Selected Matrices 里的当前矩阵
        Matrix currentMatrix;

        enum class Section { None, MatrixList, MatrixA, MatrixB };
        Section currentSection = Section::None;

        std::string line;
        while (std::getline(infile, line)) {
            line = trim(line);

            // 若空行，且在 MatrixList 区域，表示一个矩阵结束
            if (line.empty()) {
                if (currentSection == Section::MatrixList && !currentMatrix.empty()) {
                    selected_matrices.push_back(currentMatrix);
                    currentMatrix.clear();
                }
                continue;
            }

            // 检查区段切换
            if (line.find("Selected Matrices:") != std::string::npos) {
                currentSection = Section::MatrixList;
                continue;
            }
            if (line.find("Matrix A:") != std::string::npos) {
                // 若 MatrixList 部分还有未结束的矩阵，则先存入
                if (currentSection == Section::MatrixList && !currentMatrix.empty()) {
                    selected_matrices.push_back(currentMatrix);
                    currentMatrix.clear();
                }
                currentSection = Section::MatrixA;
                continue;
            }
            if (line.find("Matrix B:") != std::string::npos) {
                currentSection = Section::MatrixB;
                continue;
            }

            // 根据当前区段解析行
            std::vector<Element> rowData = parseRow(line);
            switch (currentSection) {
                case Section::MatrixList:
                    currentMatrix.push_back(rowData);
                    break;
                case Section::MatrixA:
                    matrixA.push_back(rowData);
                    break;
                case Section::MatrixB:
                    matrixB.push_back(rowData);
                    break;
                default:
                    // 未在任何区段
                    break;
            }
        }

        // 文件结束后，若处于 MatrixList 且当前矩阵未保存，则保存最后一个
        if (currentSection == Section::MatrixList && !currentMatrix.empty()) {
            selected_matrices.push_back(currentMatrix);
        }

        infile.close();
        return true;
    }
}