#include "VerilogGenerator.h"
#include "MatrixAB.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <bitset>
#include <optional>
#include <numeric>  
#include <regex>
#include "MatrixTools_for_HDL.h"

namespace MatrixProject {

std::string generate_verilog_with_matrices(const std::vector<Matrix>& matrix_list, const Matrix& matrix_a, const Matrix& matrix_b, int width, int out_width) {
    std::ostringstream verilog_code;

    int num_y = 0; // 输出个数，针对 MCM
    for (const auto& matrix : matrix_list) {
        int count = count_nonzero_diagonal(matrix);
        if (count >= num_y) {
            num_y = count;
        }
    }

    // 内部信号声明模块
    int width_l_input = 0;
    int width_r_input = 0;
    int Co_Shift_l = 0;
    int Co_Shift_r = 0;
    std::vector<int> width_v = {width};
    std::vector<int> Co_Shift_l_list;
    std::vector<int> Co_Shift_r_list;
    std::vector<int> width_l_input_list;
    std::vector<int> width_r_input_list;

    for (size_t i = 1; i < matrix_a.size(); ++i) {
        // 计算左移位量和 GCD
        std::vector<Element> Shift_l = add_column_to_list_with_list(matrix_a, i);
        // 打印 Shift__l 的内容
        print_shift(Shift_l);

        // 将 Shift__l 打平为 flat_Shift__l
        std::vector<double> flat_Shift_l = flatten_elements(Shift_l);

        if (std::all_of(flat_Shift_l.begin(), flat_Shift_l.end(), [](double s) { return std::abs(s) > 1; })) {
            Co_Shift_l = static_cast<int>(std::log2(std::abs(gcd_list(flat_Shift_l))));
        } else {
            Co_Shift_l = 0;
        }
        std::cout << "Iteration " << i << ": Co_Shift_l = " << Co_Shift_l << std::endl;

        // 计算右移位量和 GCD
        std::vector<Element> Shift_r = add_column_to_list_with_list(matrix_b, i);
        // 打印 Shift__r 的内容
        print_shift(Shift_r);

        // 将 Shift__l 打平为 flat_Shift__l
        std::vector<double> flat_Shift_r = flatten_elements(Shift_r);

        if (std::all_of(flat_Shift_r.begin(), flat_Shift_r.end(), [](double s) { return std::abs(s) > 1; })) {
            Co_Shift_r = static_cast<int>(std::log2(std::abs(gcd_list(flat_Shift_r))));
        } else {
            Co_Shift_r = 0;
        }
        std::cout << "Iteration " << i << ": Co_Shift_r = " << Co_Shift_r << std::endl;

        int max_input_width_l = max_input_width(compute_powers(Shift_l), width_v) ;
        int max_input_width_r = max_input_width(compute_powers(Shift_r), width_v) ;
        width_l_input = max_input_width_l - Co_Shift_l;
        width_r_input = max_input_width_r - Co_Shift_r;

        std::cout << "Iteration " << i << ": max_input_width_l = " << max_input_width_l << std::endl;
        std::cout << "Iteration " << i << ": max_input_width_r = " << max_input_width_r << std::endl;

        width_v.push_back(std::max(max_input_width_l, max_input_width_r) + 1);//加法后增加一位防止溢出
        Co_Shift_l_list.push_back(Co_Shift_l);
        Co_Shift_r_list.push_back(Co_Shift_r);
        width_l_input_list.push_back(width_l_input);
        width_r_input_list.push_back(width_r_input);

        verilog_code << "wire signed [" << width_v[i] -1 << ":0] v" << i << ";\n";
        verilog_code << "wire signed [" << width_l_input -1 << ":0] a_input_" << i << ";\n";
        verilog_code << "wire signed [" << width_r_input  -1<< ":0] b_input_" << i << ";\n";
    }

    /******************DEBUG*********************/
        std::cout << "width_v: [ ";
    for (const auto& val : width_v) {
        std::cout << val << " ";
    }
    std::cout << "]\n";

    std::cout << "Co_Shift_l_list: [ ";
    for (const auto& val : Co_Shift_l_list) {
        std::cout << val << " ";
    }
    std::cout << "]\n";

    std::cout << "Co_Shift_r_list: [ ";
    for (const auto& val : Co_Shift_r_list) {
        std::cout << val << " ";
    }
    std::cout << "]\n";

    std::cout << "width_l_input_list: [ ";
    for (const auto& val : width_l_input_list) {
        std::cout << val << " ";
    }
    std::cout << "]\n";

    std::cout << "width_r_input_list: [ ";
    for (const auto& val : width_r_input_list) {
        std::cout << val << " ";
    }
    std::cout << "]\n";

    //***********************中间模块*****************************
    for (size_t i = 1; i < matrix_a.size(); ++i) {
        int x = count_non_zero_elements_in_column(matrix_a, i);
        int y = count_non_zero_elements_in_column(matrix_b, i);

        // Left input
        verilog_code << "\n// v" << i << " 的左输入\n";
        if (x == 1) {
            for (size_t row = 0; row < i; ++row) {
                Element element_l = matrix_a[row][i];
                if (std::holds_alternative<double>(element_l)) {
                    double val = std::get<double>(element_l);
                    if (val != 0.0) {
                        if (val > 0) {
                            if (std::abs(val) >= 1) {
                                verilog_code << "assign a_input_" << i << " = v" << row << " <<< " << calculate_shift(val) - Co_Shift_l_list[i - 1] << ";\n";
                            } else {
                                verilog_code << "assign a_input_" << i << " = v" << row << " >>> " << calculate_shift(val) - Co_Shift_l_list[i - 1] << ";\n";
                            }
                        } else {
                            if (std::abs(val) >= 1) {
                                verilog_code << "assign a_input_" << i << " = - (v" << row << " <<< " << calculate_shift(val) - Co_Shift_l_list[i - 1] << ");\n";
                            } else {
                                verilog_code << "assign a_input_" << i << " = - (v" << row << " >>> " << calculate_shift(val) - Co_Shift_l_list[i - 1] << ");\n";
                            }
                        }
                    }
                }
            }
        } else {
            int case_num = 1;
            std::string small_expression = gen_verilog_expression(matrix_list, matrix_a, i, 0, Co_Shift_l_list);
            std::string large_expression = gen_verilog_expression(matrix_list, matrix_a, i, 1, Co_Shift_l_list);
            std::string left_expression = "sel == " + std::to_string(case_num) + " ? " + large_expression + " : " + small_expression;

            for (size_t num_matrix = 2; num_matrix < matrix_list.size(); ++num_matrix) {
                case_num++;
                std::string right_expression = gen_verilog_expression(matrix_list, matrix_a, i, num_matrix, Co_Shift_l_list);
                left_expression = "sel == " + std::to_string(case_num) + " ? " + right_expression + " : (" + left_expression + ")";
            }
            verilog_code << "assign a_input_" << i << " = " << left_expression << ";\n";
        }

        // Right input
        verilog_code << "\n// v" << i << " 的右输入\n";
        if (y == 1) {
            for (size_t row = 0; row < i; ++row) {
                Element element_r = matrix_b[row][i];
                if (std::holds_alternative<double>(element_r)) {
                    double val = std::get<double>(element_r);
                    if (val != 0.0) {
                        if (val > 0) {
                            if (std::abs(val) >= 1) {
                                verilog_code << "assign b_input_" << i << " = v" << row << " <<< " << calculate_shift(val) - Co_Shift_r_list[i - 1] << ";\n";
                            } else {
                                verilog_code << "assign b_input_" << i << " = v" << row << " >>> " << calculate_shift(val) - Co_Shift_r_list[i - 1] << ";\n";
                            }
                        } else {
                            if (std::abs(val) >= 1) {
                                verilog_code << "assign b_input_" << i << " = - (v" << row << " <<< " << calculate_shift(val) - Co_Shift_r_list[i - 1] << ");\n";
                            } else {
                                verilog_code << "assign b_input_" << i << " = - (v" << row << " >>> " << calculate_shift(val) - Co_Shift_r_list[i - 1] << ");\n";
                            }
                        }
                    }
                }
            }
        } else {
            int case_num = 1;
            std::string small_expression = gen_verilog_expression(matrix_list, matrix_b, i, 0, Co_Shift_r_list);
            std::string large_expression = gen_verilog_expression(matrix_list, matrix_b, i, 1, Co_Shift_r_list);
            std::string left_expression = "sel == " + std::to_string(case_num) + " ? " + large_expression + " : " + small_expression;

            for (size_t num_matrix = 2; num_matrix < matrix_list.size(); ++num_matrix) {
                case_num++;
                std::string right_expression = gen_verilog_expression(matrix_list, matrix_b, i, num_matrix, Co_Shift_r_list);
                left_expression = "sel == " + std::to_string(case_num) + " ? " + right_expression + " : (" + left_expression + ")";
            }
            verilog_code << "assign b_input_" << i << " = " << left_expression << ";\n";
        }

        // Compute v[i]
        verilog_code << "\n// 计算 v" << i << "\n";
        verilog_code << generate_compute_nodes(i, width_v, Co_Shift_l_list, Co_Shift_r_list, width_l_input_list, width_r_input_list);
    }

    //#*****************************输出选择部分**************************************
    int width_sel = smallest_n_for_2n_geq_x(matrix_list.size());
    verilog_code << "\n// 选择输出节点\n";

    int case_m = 0;
    std::vector<std::string> output_assignments;

    // Iterate through each matrix to generate ternary operator assignments
    for (const auto& matrix : matrix_list) {
        std::vector<int> output_index_list;

        for (size_t i = 0; i < matrix.size(); ++i) {
            if (std::holds_alternative<double>(matrix[i][i])) {
                double val = std::get<double>(matrix[i][i]);
                if (val != 0.0) {
                    std::cout << "对角线元素(单个）：";
                    print_element(matrix[i][i]);
                    output_index_list.push_back(i);
                }
            }
            else if (std::holds_alternative<std::vector<double>>(matrix[i][i])) {
                auto& vec = std::get<std::vector<double>>(matrix[i][i]);
                if (vec.size() == 1) {
                    bool has_non_zero = false;
                    for (double v : vec) {
                        if (v != 0.0) {
                            has_non_zero = true;
                            break;
                        }
                    }
                    if (has_non_zero) {
                        std::cout << "对角线元素(单个）：";
                        print_element(matrix[i][i]);
                        output_index_list.push_back(i);
                    }
                }
                else {
                    bool has_non_zero = false;
                    for (double v : vec) {
                        if (v != 0.0) {
                            has_non_zero = true;
                            break;
                        }
                    }
                    if (has_non_zero) {
                        for (int j = 0; j < vec.size(); ++j) {
                            output_index_list.push_back(i);
                        }
                        std::cout << "对角线元素（向量）：";
                        print_element(matrix[i][i]);
                        std::cout << std::endl;
                    }
                }
            }
        }

        // Pad with zeros if output node count is less than num_y
        while (output_index_list.size() < static_cast<size_t>(num_y)) {
            output_index_list.push_back(0);
        }

        /******************DEBUG*********************/
        std::cout << "output_index_list: [ ";
        for (const auto& val : output_index_list) {
            std::cout << val << " ";
        }
        std::cout << "]\n";

        // Generate ternary operator expressions for each output node
        for (int case_v = 0; case_v < num_y; ++case_v) {
            int output_index = output_index_list[case_v];
            Element out_ele = matrix[output_index][output_index];

            auto& vec = std::get<std::vector<double>>(out_ele);
            if (vec.size() == 1) {
                double val = vec[0];
                if (val != 0.0) {
                    int shift = calculate_shift(val);
                    std::string operation = (std::abs(val) >= 1) ? "<<<" : ">>>";
                    output_assignments.push_back("(sel == " + std::to_string(width_sel) + "'d" + std::to_string(case_m) + ") ? (v" + std::to_string(output_index) + " " + operation + " " + std::to_string(shift) + ")");
                } else {
                    output_assignments.push_back("(sel == " + std::to_string(width_sel) + "'d" + std::to_string(case_m) + ") ? 0");
                }
            } else {
                //auto& vec = std::get<std::vector<double>>(out_ele);
                if (vec[case_v] != 0.0) {
                    int shift = calculate_shift(vec[case_v]);
                    std::string operation = (std::abs(vec[case_v]) >= 1) ? "<<<" : ">>>";
                    output_assignments.push_back("(sel == " + std::to_string(width_sel) + "'d" + std::to_string(case_m) + ") ? (v" + std::to_string(output_index) + " " + operation + " " + std::to_string(shift) + ")");
                } else {
                    output_assignments.push_back("(sel == " + std::to_string(width_sel) + "'d" + std::to_string(case_m) + ") ? 0");
                }
            }
        }
        case_m++;
    }

    std::ostringstream header_code;
    std::vector<int> y_widths = calculate_y_width(output_assignments, width_v, num_y);

    // Generate final assignment statements for each output
    for (int case_v = 0; case_v < num_y; ++case_v) {
        std::ostringstream conditions;
        for (size_t i = 0; i < output_assignments.size(); ++i) {
            if (i % num_y == static_cast<size_t>(case_v)) {
                if (!conditions.str().empty()) {
                    conditions << " : ";
                }
                conditions << output_assignments[i];
            }
        }
        verilog_code << "assign Y" << case_v << " = " << conditions.str() << " : 0;\n";
        
    }

    // End of module
    verilog_code << "endmodule\n";

    /***********************模块定义*****************************/
    
    /******************DEBUG*********************/
    std::cout << "y_widths: [ ";
    for (const auto& val : y_widths) {
        std::cout << val << " ";
    }
    std::cout << "]\n";

    // Verilog macro definition
    header_code << "`define WIDTH " << width << "\n\n";

    // Verilog module header
    header_code << "module const_multiplier(\n";
    header_code << "  input wire signed [`WIDTH-1:0] x,\n";
    header_code << "  input wire [" << width_sel - 1 << ":0] sel,\n";

    for (int i = 0; i < num_y; ++i) {
        if (i != num_y - 1) {
            //header_code << "  output wire signed [" << y_widths[i] - 1 << ":0] y" << i << ",\n";
            header_code << "  output wire signed [" << width + out_width - 1 << ":0] Y" << i << ",\n";
        } else {
            //header_code << "  output wire signed [" << y_widths[i] - 1 << ":0] y" << i << "\n";
            header_code << "  output wire signed [" << width + out_width - 1 << ":0] Y" << i << "\n";
        }
    }
    header_code << ");\n\n";

    // Connect v0 to x
    header_code << "wire signed [`WIDTH - 1:0] v0;\n";
    header_code << "assign v0 = x ;\n";

    // Insert header_code at the beginning of verilog_code
    verilog_code.str(header_code.str() + verilog_code.str());
    return verilog_code.str();
}    
    
} // namespace MatrixProject