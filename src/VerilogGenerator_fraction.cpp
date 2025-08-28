#include "VerilogGenerator_fraction.h"
#include "MatrixAB.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <sstream>
#include <iomanip>
#include <bitset>
#include <Config.h>
#include <optional>
#include <numeric>
#include <regex>
#include "MatrixTools_for_HDL.h"
#include "MatrixCost.h"
#include "MatrixTools.h"

namespace MatrixProject {

std::string generate_verilog_with_matrices_fraction(const std::vector<Matrix>& matrix_list, const Matrix& matrix_a, const Matrix& matrix_b, int width) {
    std::ostringstream verilog_code;

    Matrix union_matrix = merge_matrices(matrix_a, matrix_b);
    print_matrix(union_matrix);

    Matrix transformed_matrix = transform_matrix(union_matrix);
    print_matrix(transformed_matrix);

    print_matrix(matrix_a);
    print_matrix(matrix_b);
    auto mac_vector = get_column_vector(transformed_matrix);
    auto max_rightshfit_vec = get_column_vector(union_matrix);

    std::vector<double> bitwidth_col;
    std::vector<double> max_rightshfit_col;

    for (size_t i = 0; i < mac_vector.size(); ++i) {
        int mac = round(abs(std::log2(mac_vector[i])));
        // An example of bit_width truncation
        if (mac > config.TRUNCATE_FACTOR){
            mac = config.TRUNCATE_FACTOR;
        }
        bitwidth_col.push_back(mac);
        double max_rightshfit = round(abs(std::log2(max_rightshfit_vec[i])));
        max_rightshfit_col.push_back(pow(2,-max_rightshfit));
    }

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
        //print_shift(Shift_l);

        // 计算右移位量和 GCD
        std::vector<Element> Shift_r = add_column_to_list_with_list(matrix_b, i);
        // 打印 Shift__r 的内容
        //print_shift(Shift_r);

        int max_input_width_l = max_input_width(compute_powers(Shift_l), width_v);
        int max_input_width_r = max_input_width(compute_powers(Shift_r), width_v);
        width_l_input = max_input_width_l;
        width_r_input = max_input_width_r;

        //std::cout << "Iteration " << i << ": max_input_width_l = " << max_input_width_l << std::endl;
        //std::cout << "Iteration " << i << ": max_input_width_r = " << max_input_width_r << std::endl;

        //width_v.push_back(std::max(max_input_width_l, max_input_width_r));
        
        width_l_input_list.push_back(width_l_input);
        width_r_input_list.push_back(width_r_input);
        
        //计算vi的位宽=每一列的最大位移+WIDTH
        //verilog_code << "wire signed [" << width_v[i] << ":0] v" << i << ";\n";//位宽增加一位防止溢出
        //width_v.push_back(width + config.TRUNCATE_FACTOR);
        width_v.push_back(width + bitwidth_col[i] );
        verilog_code << "wire signed [" << width_v[i] - 1<< ":0] v" << i << ";\n";// 在输入WIDTH的时候多加一位防止溢出

        // verilog_code << "wire signed [" << width_l_input - 1 << ":0] a_input_" << i << ";\n";
        // verilog_code << "wire signed [" << width_r_input - 1 << ":0] b_input_" << i << ";\n";
        verilog_code << "wire signed [" << width_v[i]  - 1<< ":0] a_input_" << i << ";\n"; // 在例化的时候多加一位防止溢出，其他的时不变
        verilog_code << "wire signed [" << width_v[i]  - 1<< ":0] b_input_" << i << ";\n";
    }

    /******************DEBUG*********************/
    std::cout << "bitwidth_col: [ ";
    for (const auto& val : bitwidth_col) {
        std::cout << val << " ";
    }
    std::cout << "]\n";

    // std::cout << "width_v: [ ";
    // for (const auto& val : width_v) {
    //     std::cout << val << " ";
    // }
    // std::cout << "]\n";

    // std::cout << "width_l_input_list: [ ";
    // for (const auto& val : width_v) {
    //     std::cout << val << " ";
    // }
    // std::cout << "]\n";

    // std::cout << "width_r_input_list: [ ";
    // for (const auto& val : width_v) {
    //     std::cout << val << " ";
    // }
    // std::cout << "]\n";

    //***********************中间模块*****************************
    for (size_t i = 1; i < matrix_a.size(); ++i) {
        int x = count_non_zero_elements_in_column(matrix_a, i);
        int y = count_non_zero_elements_in_column(matrix_b, i);

        // Left input
        verilog_code << "\n// v" << i << " 的左输入\n";
        if (x == 1) {
            for (size_t row = 0; row < i; ++row) {
                Element element_l = matrix_a[row][i];
                double val = std::get<double>(element_l);
                if (val != 0) {
                    std::cout << "val " << val << std::endl;
                    int right_shift = round(abs(std::log2(std::fabs(val))));
                    std::cout << "right_shift " << right_shift << std::endl;
                    int add_zeor = width_v [i] - (width_v[row] + right_shift);
                    std::cout << "add_zeor " << add_zeor << std::endl;

                    // if (val != 0.0) {
                    //     if (val > 0) {
                    //         if (std::abs(val) >= 1) {
                    //             verilog_code << "assign a_input_" << i << " = v" << row << " <<< " << calculate_shift(val)<< ";\n";
                    //         } else {
                    //             verilog_code << "assign a_input_" << i << " = v" << row << " >>> " << calculate_shift(val)<< ";\n";
                    //         }
                    //     } else {
                    //         if (std::abs(val) >= 1) {
                    //             verilog_code << "assign a_input_" << i << " = - (v" << row << " <<< " << calculate_shift(val)<< ");\n";
                    //         } else {
                    //             verilog_code << "assign a_input_" << i << " = - (v" << row << " >>> " << calculate_shift(val)<< ");\n";
                    //         }
                    //     }
                    // }

                    if (val > 0) {
                        if (add_zeor >= 0) {
                            verilog_code << "assign a_input_" << i << " = v" << row << " <<< " << add_zeor << ";\n";
                        } else {
                            verilog_code << "assign a_input_" << i << " = v" << row << " >>> " << std::abs(add_zeor) << ";\n";
                        }
                    }
                    else {
                        if (add_zeor >= 0) {
                            verilog_code << "assign a_input_" << i << " = - (v" << row << " <<< " << add_zeor << ");\n";
                        } else {
                            verilog_code << "assign a_input_" << i << " = - (v" << row << " >>> " << std::abs(add_zeor) << ");\n";
                        }
                    }

                    // if (val > 0) {
                    //     verilog_code << "assign a_input_" << i << " = v" << row << " <<< " << add_zeor << ";\n";
                    // }
                    // else {
                    //     verilog_code << "assign a_input_" << i << " = - (v" << row << " <<< " << add_zeor << ");\n";
                    // }

                }
            }
        }
        else {
            int case_num = 1;
            std::string small_expression = gen_verilog_expression_fraction(matrix_list, matrix_a, i, 0, width_v);
            std::string large_expression = gen_verilog_expression_fraction(matrix_list, matrix_a, i, 1, width_v);
            std::string left_expression = "sel == " + std::to_string(case_num) + " ? " + large_expression + " : " + small_expression;

            for (size_t num_matrix = 2; num_matrix < matrix_list.size(); ++num_matrix) {
                case_num++;
                std::string right_expression = gen_verilog_expression_fraction(matrix_list, matrix_a, i, num_matrix, width_v);
                left_expression = "sel == " + std::to_string(case_num) + " ? " + right_expression + " : (" + left_expression + ")";
            }
            verilog_code << "assign a_input_" << i << " = " << left_expression << ";\n";
        }

        // Right input
        verilog_code << "\n// v" << i << " 的右输入\n";
        if (y == 1) {
            for (size_t row = 0; row < i; ++row) {
                Element element_r = matrix_b[row][i];
                double val = std::get<double>(element_r);
                if (val != 0) {
                    std::cout << "val " << val << std::endl;
                    int right_shift = round(abs(std::log2(std::fabs(val))));
                    std::cout << "right_shift " << right_shift << std::endl;
                    int add_zeor = width_v [i] - (width_v[row] + right_shift);
                    std::cout << "add_zeor " << add_zeor << std::endl;

                    // if (val != 0.0) {
                    //     if (val > 0) {
                    //         if (std::abs(val) >= 1) {
                    //             verilog_code << "assign b_input_" << i << " = v" << row << " <<< " << calculate_shift(val)<< ";\n";
                    //         } else {
                    //             verilog_code << "assign b_input_" << i << " = v" << row << " >>> " << calculate_shift(val)<< ";\n";
                    //         }
                    //     } else {
                    //         if (std::abs(val) >= 1) {
                    //             verilog_code << "assign b_input_" << i << " = - (v" << row << " <<< " << calculate_shift(val)<< ");\n";
                    //         } else {
                    //             verilog_code << "assign b_input_" << i << " = - (v" << row << " >>> " << calculate_shift(val)<< ");\n";
                    //         }
                    //     }
                    // }

                    if (val > 0) {
                        if (add_zeor >= 0) {
                            verilog_code << "assign b_input_" << i << " = v" << row << " <<< " << add_zeor << ";\n";
                        } else {
                            verilog_code << "assign b_input_" << i << " = v" << row << " >>> " << std::abs(add_zeor) << ";\n";
                        }
                    }
                    else {
                        if (add_zeor >= 0) {
                            verilog_code << "assign b_input_" << i << " = - (v" << row << " <<< " << add_zeor << ");\n";
                        } else {
                            verilog_code << "assign b_input_" << i << " = - (v" << row << " >>> " << std::abs(add_zeor) << ");\n";
                        }
                    }
                    // if (val > 0) {
                    //     verilog_code << "assign b_input_" << i << " = v" << row << " <<< " << add_zeor << ";\n";
                    // }
                    // else {
                    //     verilog_code << "assign b_input_" << i << " = - (v" << row << " <<< " << add_zeor << ");\n";
                    // }

                }
            }
        } else {
            int case_num = 1;
            std::string small_expression = gen_verilog_expression_fraction(matrix_list, matrix_b, i, 0, width_v);
            std::string large_expression = gen_verilog_expression_fraction(matrix_list, matrix_b, i, 1, width_v);
            std::string left_expression = "sel == " + std::to_string(case_num) + " ? " + large_expression + " : " + small_expression;

            for (size_t num_matrix = 2; num_matrix < matrix_list.size(); ++num_matrix) {
                case_num++;
                std::string right_expression = gen_verilog_expression_fraction(matrix_list, matrix_b, i, num_matrix, width_v);
                left_expression = "sel == " + std::to_string(case_num) + " ? " + right_expression + " : (" + left_expression + ")";
            }
            verilog_code << "assign b_input_" << i << " = " << left_expression << ";\n";
        }

        // Compute v[i]
        verilog_code << "\n// 计算 v" << i << "\n";
        verilog_code << generate_compute_nodes_fraction(i, width_v, width_l_input_list, width_r_input_list);
    }

    //#*****************************输出选择部分**************************************
    int width_sel = smallest_n_for_2n_geq_x(matrix_list.size());
    verilog_code << "\n// 选择输出节点\n";

    int case_m = 0;
    std::vector<std::string> output_assignments;

    // Iterate through each matrix to generate ternary operator assignments
    for (const auto& matrix : matrix_list) {
        std::vector<int> output_index_list;

        // Calculate output node indices
        for (size_t i = 0; i < matrix.size(); ++i) {
            // auto& vec = std::get<std::vector<double>>(matrix[i][i]);
            if (std::holds_alternative<double>(matrix[i][i])) {
                auto& vec = std::get<double>(matrix[i][i]);
                bool has_non_zero = false;
                double val = vec;
                if (val != 0.0) {
                    has_non_zero = true;
                    break;
                }

                if (has_non_zero) {
                    std::cout << "对角线元素(单个）：";
                    print_element(matrix[i][i]);
                    output_index_list.push_back(i);
                }
            }
            // 判断元素是否为 vector<double> 并且向量中有非零元素
            else {
                auto& vec = std::get<std::vector<double>>(matrix[i][i]);
                bool has_non_zero = false;
                for (double val : vec) {
                    if (val != 0.0) {
                        has_non_zero = true;
                        break;
                    }
                }
                if (has_non_zero) {
                    for(int j=0; j < vec.size();j++){
                        output_index_list.push_back(i);
                    }
                    std::cout << "对角线元素（向量）：";
                    print_element(matrix[i][i]);
                    std::cout << std::endl;
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
                    //val = pow(2,-FRACTION_WIDTH) * pow(2,-config.TRUNCATE_FACTOR) * val;
                    // val = pow(2,-FRACTION_WIDTH) * val;
                    //int shift = calculate_shift(val);

                    int shift = width_v[output_index] - log2(val) - 1 ;
                    //std::string operation = (std::abs(val) >= 1) ? "<<<" : ">>>";
                    std::string operation;
                    if (shift > 0) {
                        operation = ">>>";
                    } else {
                        operation = "<<<";
                    }
                    output_assignments.push_back("(sel == " + std::to_string(width_sel) + "'d" + std::to_string(case_m) + ") ? (v" + std::to_string(output_index) + " " + operation + " " + std::to_string(std::abs(shift)) + ")");

                } else {
                    output_assignments.push_back("(sel == " + std::to_string(width_sel) + "'d" + std::to_string(case_m) + ") ? 0");
                }
            } else {
                //auto& vec = std::get<std::vector<double>>(out_ele);
                if (vec[case_v] != 0.0) {
                    double val = vec[case_v];
                    //val = pow(2,-FRACTION_WIDTH) * pow(2,-config.TRUNCATE_FACTOR) * vec[case_v];
                    //val = pow(2,-FRACTION_WIDTH) * vec[case_v];
                    //int shift = calculate_shift(val);
                    // int shift = width_v.back() - width_v[output_index];
                    int shift = width_v[output_index] - log2(val) - 1; //与y0对齐
                    //std::string operation = (std::abs(val) >= 1) ? "<<<" : ">>>";
                    std::string operation;

                    if (shift > 0) {
                        operation = ">>>";
                    } else {
                        operation = "<<<";
                    }
                    output_assignments.push_back("(sel == " + std::to_string(width_sel) + "'d" + std::to_string(case_m) + ") ? (v" + std::to_string(output_index) + " " + operation + " " + std::to_string(std::abs(shift)) + ")");
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
        verilog_code << "assign y" << case_v << " = " << conditions.str() << " : 0;\n";
        // if (y_widths[case_v] > width){
        //     verilog_code << "assign Y" << case_v << " = " << "y"<< case_v <<"[" <<  width - 1 << ": 0 ];\n";
        // }else{
        //     verilog_code << "assign Y" << case_v << " = " << "y"<< case_v << ";\n";
        // }
        // verilog_code << "assign Y" << case_v << " = " << "y"<< case_v <<"[" <<  width_v.back() - 1 << ": "<<  width_v.back() - 1 - width + 1 <<" ];\n";//与位宽最大的Vi对齐，截取16位
        verilog_code << "assign Y" << case_v << " = " << "y"<< case_v << ";\n";//与位宽最大的Vi对齐，截取16位

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
            header_code << "  output wire signed [" << width - 1 << ":0] Y" << i << ",\n";
        } else {
            //header_code << "  output wire signed [" << y_widths[i] - 1 << ":0] y" << i << "\n";
            header_code << "  output wire signed [" << width - 1 << ":0] Y" << i << "\n";
        }
    }
    header_code << ");\n\n";

    // Connect v0 to x
    for (int case_v = 0; case_v < num_y; ++case_v) {
        //header_code << "wire signed ["<< y_widths[case_v] <<":0] y" << case_v <<";\n";//位宽增加一位防止溢出
        header_code << "wire signed ["<< width_v.back() - 1<<":0] y" << case_v <<";\n";
    }
    //暂时取消inside-shift
    // header_code << "wire signed ["<< width + config.TRUNCATE_FACTOR <<":0] v0;\n";//位宽增加一位防止溢出
    // header_code << "assign v0 = { x[" << width - 1  << "],  x ," <<config.TRUNCATE_FACTOR<<"'b0};\n";//位宽增加一位防止溢出
    header_code << "wire signed ["<< width_v[0] - 1<<":0] v0;\n";
    header_code << "assign v0 = x ;\n";

    // Insert header_code at the beginning of verilog_code
    verilog_code.str(header_code.str() + verilog_code.str());
    return verilog_code.str();
}

} // namespace MatrixProject
