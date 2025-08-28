#ifndef VERILOG_GENERATOR_fraction_H
#define VERILOG_GENERATOR_fraction_H

#include <vector>
#include <string>
#include <optional>
#include "MatrixLoader.h"
#include "MatrixTools.h"
#include "MatrixTools_for_HDL.h"
#include <variant>

namespace MatrixProject {

std::string generate_verilog_with_matrices_fraction(const std::vector<Matrix>& matrix_list, const Matrix& matrix_a, const Matrix& matrix_b, int width);

} // namespace MatrixProject

#endif // VERILOG_GENERATOR_fraction_H
