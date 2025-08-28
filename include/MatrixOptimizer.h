#ifndef MATRIX_OPTIMIZER_H
#define MATRIX_OPTIMIZER_H

#include "MatrixLoader.h"

#include <tuple>
#include <vector>
#include <map>  // 添加此行，包含 std::map

namespace MatrixProject {

/**
 * @brief 优化矩阵组合
 * @param load_matrices 矩阵组
 * @param needPrint 是否需要打印
 * @return 选定的矩阵组、合并后的矩阵
 */
MatrixGroup optimize_matrix_combination(const MatrixGroups& load_matrices);

/**
 * @brief 动态生成遗传算法参数
 * @param load_matrices 矩阵组
 * @return 配置参数的映射
 */
std::map<std::string, double> dynamic_ga_params(const MatrixGroups& load_matrices);  // 确保函数声明以分号结尾

/**
 * @brief 遗传算法主函数
 * @param load_matrices 矩阵组
 * @param population_size 种群大小
 * @param generations 迭代代数
 * @param crossover_rate 交叉概率
 * @param mutation_rate 变异概率
 * @param elite_fraction 精英比例
 * @return 最佳解的元组，包括最佳解、最佳成本、最佳 MUX 数量
 */
std::tuple<std::vector<int>, double, double>
genetic_algorithm(const MatrixGroups& load_matrices,
                  int population_size,
                  int generations,
                  double crossover_rate,
                  double mutation_rate,
                  double elite_fraction);  // 确保函数声明以分号结尾

/**
 * @brief 用线程执行GA算法的函数
 * @param load_matrices 矩阵组
 * @return 执行结果
 */
std::vector<std::tuple<std::vector<int>, double, double>> run_ga_by_threads(const MatrixGroups& load_matrices);


}

#endif // MATRIX_OPTIMIZER_H