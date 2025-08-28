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

namespace MatrixProject {

int smallest_n_for_2n_geq_x(double x) {
    if (x <= 0.0) {
        throw std::invalid_argument("x must be greater than 0");
    }
    int n = static_cast<int>(std::ceil(std::log2(x)));
    return std::max(n, 1);
}

size_t count_non_zero_elements_in_column(const Matrix& matrix, size_t col) {
    size_t non_zero_count = 0;
    for (size_t row = 0; row < col; ++row) {
        const Element& element = matrix[row][col];

        if (std::holds_alternative<double>(element)) {
            if (std::get<double>(element) != 0.0) {
                ++non_zero_count;
            }
        } else if (std::holds_alternative<std::vector<double>>(element)) {
            const auto& vec = std::get<std::vector<double>>(element);
            for (double value : vec) {
                if (value != 0.0) {
                    ++non_zero_count;
                }
            }
        }
    }
    return non_zero_count;
}

size_t count_non_zero_elements_in_column_diag(const Matrix& matrix, size_t i) {
    size_t non_zero_count = 0;
    const Element& element = matrix[i][i];

    if (std::holds_alternative<double>(element)) {
        if (std::get<double>(element) != 0.0) {
            ++non_zero_count;
        }
    } else if (std::holds_alternative<std::vector<double>>(element)) {
        const auto& vec = std::get<std::vector<double>>(element);
        for (double value : vec) {
            if (value != 0.0) {
                ++non_zero_count;
            }
        }
    }
    return non_zero_count;
}

int calculate_shift(double value) {
    if (std::abs(value) >= 1.0) {
        return static_cast<int>(std::log2(std::abs(value)));
    } else {
        return static_cast<int>(std::log2(1.0 / std::abs(value)));
    }
}

std::vector<double> add_column_to_list(const Matrix& A, size_t i) {
    std::vector<double> list_a;
    for (size_t row = 0; row < i; ++row) {
        const Element& element = A[row][i];

        if (std::holds_alternative<double>(element)) {
            list_a.push_back(std::get<double>(element));
        } else if (std::holds_alternative<std::vector<double>>(element)) {
            const auto& vec = std::get<std::vector<double>>(element);
            list_a.insert(list_a.end(), vec.begin(), vec.end());
        }
    }
    return list_a;
}

std::vector<Element> add_column_to_list_each_row(const Matrix& A, size_t i, size_t len) {
    std::vector<Element> list_a;
    for (size_t row = 0; row < len; ++row) {
        list_a.push_back(A[row][i]);
    }
    return list_a;
}

std::vector<double> add_column_to_list_diag(const Matrix& A, size_t i) {
    std::vector<double> list_a;
    const Element& element = A[i][i];

    if (std::holds_alternative<double>(element)) {
        list_a.push_back(std::get<double>(element));
    } else if (std::holds_alternative<std::vector<double>>(element)) {
        const auto& vec = std::get<std::vector<double>>(element);
        list_a.insert(list_a.end(), vec.begin(), vec.end());
    }
    return list_a;
}

size_t count_elements_before_index(const std::vector<Element>& lst, size_t i) {
    size_t count = 0;
    for (size_t j = 0; j < i; ++j) {
        const Element& elem = lst[j];
        if (std::holds_alternative<double>(elem)) {
            ++count;
        } else if (std::holds_alternative<std::vector<double>>(elem)) {
            count += std::get<std::vector<double>>(elem).size();
        }
    }
    return count;
}

std::optional<size_t> find_element_in_list_each_row(
    const Matrix& matrix,
    const std::vector<double>& list_vals,
    const std::vector<Element>& list_each_row,
    size_t col,
    size_t row
) {
    const Element& element = matrix[row][col];

    if (std::holds_alternative<double>(element)) {
        double elem_value = std::get<double>(element);

        if (std::holds_alternative<double>(list_each_row[row])) {
            double row_value = std::get<double>(list_each_row[row]);
            if (elem_value == row_value) {
                return count_elements_before_index(list_each_row, row);
            }
        } else if (std::holds_alternative<std::vector<double>>(list_each_row[row])) {
            const auto& vec = std::get<std::vector<double>>(list_each_row[row]);
            auto it = std::find(vec.begin(), vec.end(), elem_value);
            if (it != vec.end()) {
                size_t index_in_row = std::distance(vec.begin(), it);
                return index_in_row + count_elements_before_index(list_each_row, row);
            }
        }
    } else if (std::holds_alternative<std::vector<double>>(element)) {
        const auto& elem_vec = std::get<std::vector<double>>(element);
        for (double elem_value : elem_vec) {
            if (std::holds_alternative<double>(list_each_row[row])) {
                double row_value = std::get<double>(list_each_row[row]);
                if (elem_value == row_value) {
                    return count_elements_before_index(list_each_row, row);
                }
            } else if (std::holds_alternative<std::vector<double>>(list_each_row[row])) {
                const auto& vec = std::get<std::vector<double>>(list_each_row[row]);
                auto it = std::find(vec.begin(), vec.end(), elem_value);
                if (it != vec.end()) {
                    size_t index_in_row = std::distance(vec.begin(), it);
                    return index_in_row + count_elements_before_index(list_each_row, row);
                }
            }
        }
    }
    return std::nullopt;
}

std::string int_to_binary_fixed_width(int n, int width) {
    std::string binary_str = std::bitset<64>(n).to_string();
    // 取右边的 width 位
    return binary_str.substr(64 - width, width);
}

std::vector<std::string> pad_binary_list(const std::vector<std::string>& binary_list) {
    // 找到最大位宽
    size_t max_width = 0;
    for (const auto& bin_str : binary_list) {
        if (bin_str.length() > max_width) {
            max_width = bin_str.length();
        }
    }

    // 对位宽较小的元素进行补0操作
    std::vector<std::string> padded_list;
    for (const auto& bin_str : binary_list) {
        std::string padded_str = bin_str;
        padded_str.append(max_width - bin_str.length(), '0');
        padded_list.push_back(padded_str);
    }

    return padded_list;
}

int max_nonzero_diagonal_index(const Matrix& matrix) {
    int n = static_cast<int>(matrix.size());
    int max_index = -1;

    for (int i = 0; i < n; ++i) {
        const Element& elem = matrix[i][i];
        bool is_nonzero = false;

        if (std::holds_alternative<double>(elem)) {
            is_nonzero = std::get<double>(elem) != 0.0;
        } else if (std::holds_alternative<std::vector<double>>(elem)) {
            const auto& vec = std::get<std::vector<double>>(elem);
            is_nonzero = std::any_of(vec.begin(), vec.end(), [](double v) { return v != 0.0; });
        }

        if (is_nonzero) {
            max_index = i;
        }
    }

    return max_index;
}

size_t count_nonzero_diagonal(const Matrix& matrix) {
    size_t count = 0;
    size_t n = matrix.size();

    for (size_t i = 0; i < n; ++i) {
        const Element& elem = matrix[i][i];

        if (std::holds_alternative<double>(elem)) {
            if (std::get<double>(elem) != 0.0) {
                ++count;
            }
        } else if (std::holds_alternative<std::vector<double>>(elem)) {
            const auto& vec = std::get<std::vector<double>>(elem);
            if (std::any_of(vec.begin(), vec.end(), [](double v) { return v != 0.0; })) {
                ++count;
            }
        }
    }

    return count;
}

// 
std::vector<Element> add_column_to_list_with_list(const Matrix& A, size_t i) {
    std::vector<Element> list_a;
    for (size_t row = 0; row < i; ++row) {
        Element element = A[row][i];

        // 检查元素是否为非零
        list_a.push_back(element);
        
    }
    return list_a;
}

std::vector<std::variant<int, std::vector<int>>> compute_powers(const std::vector<Element>& shift_list) {
    auto get_power_of_two = [](const Element& value) -> std::variant<int, std::vector<int>> {
        if (std::holds_alternative<double>(value)) {
            double val = std::get<double>(value);
            return (val == 0) ? 0 : static_cast<int>(std::log2(std::abs(val)));
        } else if (std::holds_alternative<std::vector<double>>(value)) {
            const auto& vec = std::get<std::vector<double>>(value);
            std::vector<int> result;
            for (double v : vec) {
                result.push_back((v == 0) ? 0 : static_cast<int>(std::log2(std::abs(v))));
            }
            return result;
        } else {
            throw std::invalid_argument("Unsupported data type in shift_list. Only doubles and vectors of doubles are allowed.");
        }
    };

    std::vector<std::variant<int, std::vector<int>>> powers;
    for (const auto& item : shift_list) {
        powers.push_back(get_power_of_two(item));
    }

    return powers;
}

int max_input_width(const std::vector<std::variant<int, std::vector<int>>>& new_shift, const std::vector<int>& width_v) {
    std::vector<int> result;

    for (size_t i = 0; i < new_shift.size(); ++i) {
        if (std::holds_alternative<int>(new_shift[i])) {
            int ns = std::get<int>(new_shift[i]);
            result.push_back(ns + width_v[i]);
        } else if (std::holds_alternative<std::vector<int>>(new_shift[i])) {
            const auto& ns_list = std::get<std::vector<int>>(new_shift[i]);
            int max_val = *std::max_element(ns_list.begin(), ns_list.end());
            result.push_back(width_v[i] + max_val);
        } else {
            throw std::invalid_argument("Unsupported data type in new_shift. Only integers and vectors of integers are allowed.");
        }
    }

    return *std::max_element(result.begin(), result.end());
}

int gcd_list(const std::vector<double>& list) {
    // Convert floating-point numbers to integers by rounding
    std::vector<int> int_lst;
    for (double s : list) {
        int_lst.push_back(static_cast<int>(abs(std::round(s))));
    }

    // Calculate the GCD of the list
    return std::reduce(int_lst.begin(), int_lst.end(), int_lst[0], std::gcd<int, int>);
}

std::string gen_verilog_expression(const MatrixGroup& matrix_list, const Matrix& matrix_a, size_t i, size_t num_matrix, const std::vector<int>& Co_Shift_list) {
    for (size_t row = 0; row < i; ++row) {
        Element element_l = matrix_a[row][i];
        if (num_matrix >= matrix_list.size()) {
            return " ";
        }
        Element corresponding_element = matrix_list[num_matrix][row][i];

        // Ensure element_l and corresponding_element are lists
        std::vector<double> element_l_list;
        if (std::holds_alternative<double>(element_l)) {
            element_l_list.push_back(std::get<double>(element_l));
        } else {
            element_l_list = std::get<std::vector<double>>(element_l);
        }

        std::vector<double> corresponding_element_list;
        if (std::holds_alternative<double>(corresponding_element)) {
            corresponding_element_list.push_back(std::get<double>(corresponding_element));
        } else {
            corresponding_element_list = std::get<std::vector<double>>(corresponding_element);
        }

        // Check for matching elements
        for (double el : element_l_list) {
            for (double ce : corresponding_element_list) {
                if (el == ce && std::abs(el) > 0) {
                    if (el > 0) {
                        return (std::abs(el) >= 1) ? 
                            "v" + std::to_string(row) + " <<< " + std::to_string(calculate_shift(el) - Co_Shift_list[i - 1]) :
                            "v" + std::to_string(row) + " >>> " + std::to_string(calculate_shift(el) - Co_Shift_list[i - 1]);
                    } else {
                        return (std::abs(el) >= 1) ? 
                            "- (v" + std::to_string(row) + " <<< " + std::to_string(calculate_shift(el) - Co_Shift_list[i - 1]) + ")" :
                            "- (v" + std::to_string(row) + " >>> " + std::to_string(calculate_shift(el) - Co_Shift_list[i - 1]) + ")";
                    }
                }
            }
        }
    }

    return " 0 ";
}

std::string gen_verilog_expression_fraction(const MatrixGroup& matrix_list, const Matrix& matrix_a, size_t i, size_t num_matrix, std::vector<int> width_v) {
    for (size_t row = 0; row < i; ++row) {
        Element element_l = matrix_a[row][i];
        if (num_matrix >= matrix_list.size()) {
            return " ";
        }
        Element corresponding_element = matrix_list[num_matrix][row][i];

        // Ensure element_l and corresponding_element are lists
        std::vector<double> element_l_list;
        if (std::holds_alternative<double>(element_l)) {
            element_l_list.push_back(std::get<double>(element_l));
        } else {
            element_l_list = std::get<std::vector<double>>(element_l);
        }

        std::vector<double> corresponding_element_list;
        if (std::holds_alternative<double>(corresponding_element)) {
            corresponding_element_list.push_back(std::get<double>(corresponding_element));
        } else {
            corresponding_element_list = std::get<std::vector<double>>(corresponding_element);
        }

        // Check for matching elements
        for (double el : element_l_list) {
            for (double ce : corresponding_element_list) {
                if (el == ce && std::abs(el) > 0) {
                        // if (el > 0) {
                        //     if (std::abs(el) >= 1) {
                        //             //verilog_code << "assign a_input_" << i << " = v" << row << " <<< " << calculate_shift(val)<< ";\n";
                        //             return  "v" + std::to_string(row) + " <<< " + std::to_string(calculate_shift(el));
                        //         } else {
                        //             //verilog_code << "assign a_input_" << i << " = v" << row << " >>> " << calculate_shift(val)<< ";\n";
                        //             return  "v" + std::to_string(row) + " >>> " + std::to_string(calculate_shift(el));
                        //         }
                        //     //return  "v" + std::to_string(row) + " >>> " + std::to_string(calculate_shift(el));
                        // } else {
                        //     if (std::abs(el) >= 1) {
                        //             //verilog_code << "assign a_input_" << i << " = v" << row << " <<< " << calculate_shift(val)<< ";\n";
                        //             return  "- (v" + std::to_string(row) + " <<< " + std::to_string(calculate_shift(el)) + ")";
                        //         } else {
                        //             //verilog_code << "assign a_input_" << i << " = v" << row << " >>> " << calculate_shift(val)<< ";\n";
                        //             return  "- (v" + std::to_string(row) + " >>> " + std::to_string(calculate_shift(el)) + ")";
                        //         }
                        //     //return   "- (v" + std::to_string(row) + " >>> " + std::to_string(calculate_shift(el)) + ")";
                        // }

                        std::cout << "val " << el << std::endl;
                        int right_shift = round(abs(std::log2(std::fabs(el))));
                        std::cout << "right_shift " << right_shift << std::endl;
                        int add_zeor = width_v [i] - (width_v[row] + right_shift);
                        std::cout << "add_zeor " << add_zeor << std::endl;
                        // if (el > 0) {
                        //     return  "v" + std::to_string(row) + " <<< " + std::to_string(add_zeor);
                        // }
                        // else {
                        //     return  "- (v" + std::to_string(row) + " <<< " + std::to_string(add_zeor) + ")";
                        // }
                        if (el > 0) {
                            if (add_zeor >= 0) {
                                return  "v" + std::to_string(row) + " <<< " + std::to_string(add_zeor);
                            } else {
                                return  "v" + std::to_string(row) + " >>> " + std::to_string(std::abs(add_zeor));
                            }
                        }
                        else {
                            if (add_zeor >= 0) {
                                return  "- (v" + std::to_string(row) + " <<< " + std::to_string(add_zeor) + ")";
                            } else {
                                return  "- (v" + std::to_string(row) + " >>> " + std::to_string(std::abs(add_zeor)) + ")";
                            }
                        }
                    }
                }
            }
        
    }

    return " 0 ";
}

std::string generate_compute_nodes(
    int i, 
    const std::vector<int>& width_v, 
    const std::vector<int>& Co_Shift_l, 
    const std::vector<int>& Co_Shift_r,
    const std::vector<int>& width_l_input_list,
    const std::vector<int>& width_r_input_list
) {
    std::ostringstream verilog_code;
    std::string l_input = "a_input_" + std::to_string(i);
    std::string r_input = "b_input_" + std::to_string(i);
    std::ostringstream vi;

    if (Co_Shift_l[i-1] == 0 && Co_Shift_r[i-1] == 0){//当共同位移都为0
        int diff_w_l = width_v[i] - width_l_input_list[i-1];
        int diff_w_r = width_v[i] - width_r_input_list[i-1];
        std::string l_input_expression = "{{{" + std::to_string(diff_w_l) + "{a_input_" + std::to_string(i) + "[" + std::to_string((width_l_input_list[i-1] - 1)) + "]}}},"+l_input+ "}";
        std::string r_input_expression = "{{{" + std::to_string(diff_w_r) + "{b_input_" + std::to_string(i) + "[" + std::to_string((width_r_input_list[i-1] - 1)) + "]}}},"+r_input+ "};";
        vi << "assign v" << i << "[" << (width_v[i] - 1) << ": 0 ] = " << l_input_expression << "+" <<r_input_expression;
    }else if (Co_Shift_l[i-1] >= Co_Shift_r[i-1]) {//当左输入位移多一些
        if (Co_Shift_r[i-1] > 0) {//当右输入位移不为0
            vi << "assign v" << i << "[" << (Co_Shift_r[i-1] - 1) << ":0] = 0;\n";
            if (Co_Shift_l[i-1] > Co_Shift_r[i-1]){//当右输入位移==左位移时，跳过
                vi << "assign v" << i << "[" << (Co_Shift_l[i-1] - 1) << ":" << Co_Shift_r[i-1] << "] = " 
               << r_input << "[" << (Co_Shift_l[i-1] - Co_Shift_r[i-1] - 1) << ":0];\n";
            }

            int diff_w_l = width_v[i] - (width_l_input_list[i-1] + Co_Shift_l[i-1]);
            int diff_w_r = width_v[i] - (width_r_input_list[i-1] + Co_Shift_r[i-1]);
            std::string l_input_expression = "{{{" + std::to_string(diff_w_l) + "{a_input_" + std::to_string(i) + "[" + std::to_string((width_l_input_list[i-1] - 1)) + "]}}},"+l_input+ "}";
            std::string r_input_expression = "{{{" + std::to_string(diff_w_r) + "{b_input_" + std::to_string(i) + "[" + std::to_string((width_r_input_list[i-1] - 1)) + "]}}},"+r_input+"[" + std::to_string((width_r_input_list[i-1] - 1)) +":"+ std::to_string((Co_Shift_l[i-1] - Co_Shift_r[i-1]))+ "]};";
            vi << "assign v" << i << "[" << (width_v[i] - 1) << ":" << Co_Shift_l[i-1] << "] = " << l_input_expression << "+" <<r_input_expression;
        } else {//当右输入位移为0
            vi << "assign v" << i << "[" << (Co_Shift_l[i-1] - 1) << ":" << Co_Shift_r[i-1] << "] = " 
               << r_input << "[" << (Co_Shift_l[i-1] - Co_Shift_r[i-1] - 1) << ":0];\n";

            int diff_w_l = width_v[i] - (width_l_input_list[i-1] + Co_Shift_l[i-1]);
            int diff_w_r = width_v[i] - (width_r_input_list[i-1] + Co_Shift_r[i-1]);
            std::string l_input_expression = "{{{" + std::to_string(diff_w_l) + "{a_input_" + std::to_string(i) + "[" + std::to_string((width_l_input_list[i-1] - 1)) + "]}}},"+l_input+ "}";
            std::string r_input_expression = "{{{" + std::to_string(diff_w_r) + "{b_input_" + std::to_string(i) + "[" + std::to_string((width_r_input_list[i-1] - 1)) + "]}}},"+r_input+"[" + std::to_string((width_r_input_list[i-1] - 1)) +":"+ std::to_string((Co_Shift_l[i-1] - Co_Shift_r[i-1]))+ "]};";
            vi << "assign v" << i << "[" << (width_v[i] - 1) << ":" << Co_Shift_l[i-1] << "] = " << l_input_expression << "+" <<r_input_expression;
        }
    } else {//当右输入位移多一些
        if (Co_Shift_l[i-1] > 0) {//当左输入位移不为0
            vi << "assign v" << i << "[" << (Co_Shift_l[i-1] - 1) << ":0] = 0;\n";
            if (Co_Shift_r[i-1] > Co_Shift_l[i-1]){//当右输入位移==左位移时，跳过
                vi << "assign v" << i << "[" << (Co_Shift_r[i-1] - 1) << ":" << Co_Shift_l[i-1] << "] = " 
                << l_input << "[" << (Co_Shift_r[i-1] - Co_Shift_l[i-1] - 1) << ":0];\n";
            }

            int diff_w_l = width_v[i] - (width_l_input_list[i-1] + Co_Shift_l[i-1]);
            int diff_w_r = width_v[i] - (width_r_input_list[i-1] + Co_Shift_r[i-1]);
            std::string l_input_expression = "{{{" + std::to_string(diff_w_l) + "{a_input_" + std::to_string(i) + "[" + std::to_string((width_l_input_list[i-1] - 1)) + "]}}},"+l_input+"[" + std::to_string((width_l_input_list[i-1] - 1)) +":"+ std::to_string((Co_Shift_r[i-1] - Co_Shift_l[i-1]))+ "]}";
            std::string r_input_expression = "{{{" + std::to_string(diff_w_r) + "{b_input_" + std::to_string(i) + "[" + std::to_string((width_r_input_list[i-1] - 1)) + "]}}},"+r_input+ "};";
            vi << "assign v" << i << "[" << (width_v[i] - 1) << ":" << Co_Shift_r[i-1] << "] = " << l_input_expression << "+" <<r_input_expression;
        } else {//当左输入位移为0
            vi << "assign v" << i << "[" << (Co_Shift_r[i-1] - 1) << ":" << Co_Shift_l[i-1] << "] = " 
               << l_input << "[" << (Co_Shift_r[i-1] - Co_Shift_l[i-1] - 1) << ":0];\n";
            
            int diff_w_l = width_v[i] - (width_l_input_list[i-1] + Co_Shift_l[i-1]);
            int diff_w_r = width_v[i] - (width_r_input_list[i-1] + Co_Shift_r[i-1]);
            std::string l_input_expression = "{{{" + std::to_string(diff_w_l) + "{a_input_" + std::to_string(i) + "[" + std::to_string((width_l_input_list[i-1] - 1)) + "]}}},"+l_input+"[" + std::to_string((width_l_input_list[i-1] - 1)) +":"+ std::to_string(Co_Shift_r[i-1] - Co_Shift_l[i-1])+ "]}";
            std::string r_input_expression = "{{{" + std::to_string(diff_w_r) + "{b_input_" + std::to_string(i) + "[" + std::to_string((width_r_input_list[i-1] - 1)) + "]}}},"+r_input+ "};";
            vi << "assign v" << i << "[" << (width_v[i] - 1) << ":" << Co_Shift_r[i-1] << "] = " << l_input_expression << "+" <<r_input_expression;
        }
    }
    verilog_code << vi.str();
    return verilog_code.str();
}

std::string generate_compute_nodes_fraction(
    int i, 
    const std::vector<int>& width_v, 
    const std::vector<int>& width_l_input_list,
    const std::vector<int>& width_r_input_list
) {
    std::ostringstream verilog_code;
    std::string l_input = "a_input_" + std::to_string(i);
    std::string r_input = "b_input_" + std::to_string(i);
    std::ostringstream vi;

    if ( width_l_input_list[i-1] > width_r_input_list[i-1] ) {
        int diff_w = width_l_input_list[i-1] - width_r_input_list[i-1];
        vi << "assign v" << i << " = " << l_input 
            // << " + {{{" << diff_w << "{" << "b_input_" << i << "[" << (width_r_input_list[i-1] - 1) << "]}}}," 
            // << r_input << "};\n";
            <<  " + " << r_input << ";\n";
    } else {
        int diff_w = width_r_input_list[i-1] - width_l_input_list[i-1];
        vi << "assign v" << i <<" = " << r_input 
            // << " + {{{" << diff_w << "{" << "a_input_" << i << "[" << (width_l_input_list[i-1] - 1) << "]}}}," 
            // << l_input <<"};\n";
            <<  " + " <<l_input <<";\n";
    }
    verilog_code << vi.str();
    return verilog_code.str();
}

std::vector<int> calculate_y_width(const std::vector<std::string>& output_assignments, const std::vector<int>& v_widths, int num_y) {
    std::vector<int> y_widths;

    // Regular expression to match v{index} and shift operations
    std::regex pattern(R"(v(\d+)\s*([<>]{3})\s*(\d+))");

    for (int case_v = 0; case_v < num_y; ++case_v) {
        int max_width = 0;

        // Extract expressions related to y{case_v}
        std::vector<std::string> related_exprs;
        for (size_t i = 0; i < output_assignments.size(); ++i) {
            if (i % num_y == case_v) {
                related_exprs.push_back(output_assignments[i]);
            }
        }

        for (const auto& expr : related_exprs) {
            std::smatch match;
            if (std::regex_search(expr, match, pattern)) {
                int v_index = std::stoi(match[1].str());
                std::string shift_op = match[2].str();
                int shift_value = std::stoi(match[3].str());

                // Get the width of the v signal
                if (v_index < v_widths.size()) {
                    int v_width = v_widths[v_index];

                    // Calculate the resulting width based on the shift operation
                    int result_width = (shift_op == "<<<") ? v_width + shift_value : v_width - shift_value ;

                    max_width = std::max(max_width, result_width);
                } else {
                    throw std::out_of_range("v" + std::to_string(v_index) + " index out of range in v_widths");
                }
            }
        }

        y_widths.push_back(max_width);
    }

    return y_widths;
}

// 打印 Element 类型的辅助函数
void print_element(const Element& elem) {
    if (std::holds_alternative<double>(elem)) {
        double val = std::get<double>(elem);
        std::cout << ((val == 0) ? 0 : static_cast<int>(std::log2(std::abs(val)))) << " ";
    } else if (std::holds_alternative<std::vector<double>>(elem)) {
        std::cout << "[";
        for (const auto& val : std::get<std::vector<double>>(elem)) {
            std::cout << ((val == 0) ? 0 : static_cast<int>(std::log2(std::abs(val)))) << " ";
        }
        std::cout << "] ";
    } else {
        std::cout << "Unknown Element ";
    }
}

// 打印 std::vector<Element> 类型的辅助函数
void print_shift(const std::vector<Element>& shift_l) {
    std::cout << "Shift = [";
    for (const auto& elem : shift_l) {
        print_element(elem);
    }
    std::cout << "]" << std::endl;
}

// Flatten Elements and keep only non-zero values
std::vector<double> flatten_elements(const std::vector<Element>& elements) {
    std::vector<double> flat_list;
    for (const auto& elem : elements) {
        if (std::holds_alternative<double>(elem)) {
            double val = std::get<double>(elem);
            if (val != 0.0) {
                flat_list.push_back(val);
            }
        } else if (std::holds_alternative<std::vector<double>>(elem)) {
            const auto& vec = std::get<std::vector<double>>(elem);
            for (double val : vec) {
                if (val != 0.0) {
                    flat_list.push_back(val);
                }
            }
        }
    }
    return flat_list;
}
 
 
} // namespace MatrixProject