#include "MatrixOptimizer.h"
#include "MatrixTools.h"
#include "MatrixCost.h"
#include "MatrixAB.h"
#include "thread.h"

#include <random>
#include <algorithm>
#include <tuple>
#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <thread>
#include "Config.h"
#include <fstream>
#include <ctime>
#include <filesystem>
#include <sstream>  // 用于 std::ostringstream
#include <queue>
#include <unordered_set>


namespace MatrixProject {

// -------------------- Top‑K 穷举辅助结构 --------------------
// 说明：本节代码无需依赖遗传算法即可从所有组合中找出
//       Jaccard 相似度最高的前 k 个解。

// ResultIdx：存储每个 MatrixGroup 选中矩阵的索引路径
using ResultIdx = std::vector<int>;

// ResultEntry：用于在优先队列中保存一个候选解
//   score    —— 该组合的全局 Jaccard 相似度
//   indices  —— 选择路径（与 groups 的顺序一一对应）
//   mux_num  —— 预留字段，保持与 GA 版本返回值一致
struct ResultEntry {
    double     score;      // Jaccard 分数
    ResultIdx  indices;    // 组合的索引向量
    double     mux_num = 0.0;
};

// MinCmp：最小堆比较器，使得优先队列 top 顶端始终为
//         “当前第 k 名” 的最小分数，可在插入时自动淘汰
struct MinCmp {
    bool operator()(const ResultEntry& a, const ResultEntry& b) const {
        return a.score > b.score;   // 分数小的优先弹出
    }
};

// 将 Element 转为集合，并忽略所有值为 0 的元素（忽略正负号，取绝对值）
std::unordered_set<double> elem_to_set(const Element& e)
{
    std::unordered_set<double> s;

    std::visit([&](auto&& arg) {
        using T = std::decay_t<decltype(arg)>;
        if constexpr (std::is_same_v<T, double>) {
            s.insert(std::fabs(arg));                // 单个 double，取绝对值
        } else {
            for (double v : arg)                     // vector<double>
                s.insert(std::fabs(v));              // 逐个取绝对值
        }
    }, e);

    s.insert(0.0);   // 规则：所有位置都默认含有 0
    return s;        // 若原来已有 0 不受影响；若没有就新增
}

// ------------------------------------------------------------

// 计算若干矩阵的“全局 Jaccard 相似度”
// 步骤：遍历每个坐标位置 -> 求所有矩阵在该坐标的集合交/并集 ->
//       累加交集大小与并集大小，最后取比值
double global_jaccard(const std::vector<const Matrix*>& mats)
{
    std::size_t rows = mats[0]->size();
    std::size_t cols = (*mats[0])[0].size();

    double sum_J = 0.0;                 // 累加所有格子的 Jaccard

    for (std::size_t i = 0; i < rows; ++i)
        for (std::size_t j = 0; j < cols; ++j)
        {
            // ① 先取第 1 张矩阵的集合作初始交/并集
            std::unordered_set<double> inter = elem_to_set((*mats[0])[i][j]);
            std::unordered_set<double> uni   = inter;   // 拷贝一份

            // ② 与后续矩阵在同坐标的集合做交并
            for (std::size_t k = 1; k < mats.size(); ++k) {
                auto s = elem_to_set((*mats[k])[i][j]);

                // 并集更新
                uni.insert(s.begin(), s.end());

                // 交集更新：留下同时存在的元素
                std::unordered_set<double> tmp;
                for (double v : inter)
                    if (s.count(v)) tmp.insert(v);
                inter.swap(tmp);
            }

            // ③ 单格 Jaccard
            double Jij = uni.empty() ? 1.0                       // 全为空集合时视为完全一致
                                     : static_cast<double>(inter.size()) /
                                       static_cast<double>(uni.size());

            sum_J += Jij;        // ④ 累加
        }

    return sum_J;    // 返回所有坐标 Jaccard 之和（范围 0 ~ rows*cols）
}

// 深度优先枚举每个 MatrixGroup 的选择组合，并实时维护 Top‑K
void dfs_enum(const MatrixGroups& groups, std::size_t depth, ResultIdx& idx_buf, std::vector<const Matrix*>& mats_buf,
              std::priority_queue<ResultEntry, std::vector<ResultEntry>, MinCmp>& top, int k)
{
    if (depth == groups.size()) {
        // 已经选完一条路径，计算分数并尝试放入堆
        double score = global_jaccard(mats_buf);

        // ----------- 调试输出：打印当前组合及其 Jaccard 分数 -----------
        MatrixGroup current_group;
        for (const Matrix* mp : mats_buf) current_group.push_back(*mp);

        std::cout << "=== Combination evaluated, Jaccard = " << score << " ===\n";
        print_matrix(current_group);               // 已有的 MatrixGroup 打印函数
        std::cout << "---------------------------------------------\n";

        if (top.size() < static_cast<std::size_t>(k))
            top.push({score, idx_buf});
        else if (score > top.top().score) {  // 淘汰当前最差
            top.pop();
            top.push({score, idx_buf});
        }
        return;
    }

    // 递归枚举 depth 这一层的所有候选矩阵
    const auto& g = groups[depth];
    for (std::size_t i = 0; i < g.size(); ++i) {
        idx_buf[depth]  = static_cast<int>(i);
        mats_buf[depth] = &g[i];
        dfs_enum(groups, depth + 1, idx_buf, mats_buf, top, k);
    }
}

// 对外接口：返回与 GA 兼容的 (indices, score, mux_num) 向量
std::vector<std::tuple<std::vector<int>, double, double>> exhaustive_top_k(const MatrixGroups& groups, int k /*默认取 10 */)
{
    if (groups.empty()) return {};

    // Top‑K 使用最小堆维护
    std::priority_queue<ResultEntry, std::vector<ResultEntry>, MinCmp> top;

    ResultIdx idx_buf(groups.size());                  // 当前路径
    std::vector<const Matrix*> mats_buf(groups.size()); // 当前矩阵指针

    dfs_enum(groups, 0, idx_buf, mats_buf, top, k);

    // 将堆中结果倒序输出（得分从高到低）
    std::vector<std::tuple<std::vector<int>, double, double>> out;
    while (!top.empty()) {
        const auto& e = top.top();
        out.emplace_back(e.indices, e.score, e.mux_num);
        top.pop();
    }
    std::reverse(out.begin(), out.end());
    return out;
}

std::map<std::string, double> dynamic_ga_params(const MatrixGroups& load_matrices) {
    // 计算解空间大小
    double total_combinations = 1;
    for (const auto& group : load_matrices) {
        total_combinations *= group.size();
    }

    // 使用对数缩放增强算法参数
    int population_size = std::min(4000, std::max(1000, static_cast<int>(300 * std::log10(total_combinations / 10'000'000 + 10))));
    int generations = std::min(6000, std::max(1500, static_cast<int>(700 * std::log10(total_combinations / 10'000'000 + 10))));
    double crossover_rate = std::max(0.85, 0.95 - 0.03 * std::log10(total_combinations / 1'000'000'000 + 1));
    double max_mutation_rate = std::min(0.6, std::max(0.3, 0.3 + 0.15 * std::log10(total_combinations / 1'000'000'000 + 1)));
    int tournament_size = std::max(15, static_cast<int>(0.07 * population_size));
    double elite_fraction = std::max(0.02, 0.06 - 0.02 * std::log10(total_combinations / 1'000'000'000 + 1));

    std::map<std::string, double> GAConfig = {
        {"population_size", static_cast<double>(population_size)},
        {"generations", static_cast<double>(generations)},
        {"crossover_rate", crossover_rate},
        {"max_mutation_rate", max_mutation_rate},
        {"tournament_size", static_cast<double>(tournament_size)},
        {"elite_fraction", elite_fraction}
    };

    std::cout << "Enhanced Genetic Algorithm Parameters for Large Scale Optimization:" << std::endl;
    for (const auto& [key, value] : GAConfig) {
        std::cout << key << ": " << value << std::endl;
    }

    return GAConfig;
}

std::tuple<std::vector<int>, double, double>
genetic_algorithm(const MatrixGroups& load_matrices,
                  int population_size,
                  int generations,
                  double crossover_rate,
                  double mutation_rate,
                  double elite_fraction) {
    // 计算精英保留的个体数量
    int elite_size = std::max(1, static_cast<int>(population_size * elite_fraction));

    // 初始化随机数生成器
    std::random_device rd;
    std::mt19937 gen(rd());

    // 初始化种群
    std::vector<std::vector<int>> population(population_size);
    for (auto& individual : population) {
        individual.resize(load_matrices.size());
        for (size_t i = 0; i < load_matrices.size(); ++i) {
            std::uniform_int_distribution<> dis(0, static_cast<int>(load_matrices[i].size()) - 1);
            individual[i] = dis(gen);
        }
    }

    // 定义评估函数
    auto evaluate = [&](const std::vector<int>& individual) -> std::tuple<double, double, double> {
        // 计算个体适应度
        MatrixGroup selected_matrices;
        for (size_t i = 0; i < individual.size(); ++i) {
            selected_matrices.push_back(load_matrices[i][individual[i]]);
        }

        // 使用 merge_matrices_custom 函数合并矩阵
        Matrix A, B;
        std::tie(A, B) = merge_matrices_custom(selected_matrices, false);

        double cost = 0.0;
        double MUX_num = 0.0;
        std::tie(cost, MUX_num) = cost_analysis(A, B);

        // 计算加权适应度值
        double fitness = 0.7 * cost + 0.3 * MUX_num;
        return std::make_tuple(cost, MUX_num, fitness);
    };


    // 定义选择父代的函数（锦标赛选择法）
    auto select_parents = [&](const std::vector<std::vector<int>>& population, const std::vector<std::tuple<double, double, double>>& scores) -> std::vector<std::vector<int>> {
        int tournament_size = 3;
        std::vector<std::vector<int>> selected;
        std::uniform_int_distribution<> dis(0, population_size - 1);

        for (int i = 0; i < 2; ++i) {
            std::vector<std::pair<std::vector<int>, double>> tournament;
            for (int j = 0; j < tournament_size; ++j) {
                int idx = dis(gen);
                tournament.emplace_back(population[idx], std::get<2>(scores[idx])); // 使用加权适应度值进行比较
            }
            auto winner = *std::min_element(tournament.begin(), tournament.end(), [](const auto& a, const auto& b) {
                return a.second < b.second;
            });
            selected.push_back(winner.first);
        }
        return selected;
    };

    // 定义交叉操作
    auto crossover = [&](const std::vector<int>& parent1, const std::vector<int>& parent2) -> std::vector<int> {
        std::uniform_real_distribution<> dis(0.0, 1.0);
        if (dis(gen) < crossover_rate) {
            std::uniform_int_distribution<> dis_point(0, static_cast<int>(parent1.size()) - 1);
            int crossover_point = dis_point(gen);
            std::vector<int> child(parent1.begin(), parent1.begin() + crossover_point);
            child.insert(child.end(), parent2.begin() + crossover_point, parent2.end());
            return child;
        }
        return dis(gen) < 0.5 ? parent1 : parent2;
    };

    // 定义变异操作
    auto mutate = [&](std::vector<int>& individual) {
        std::uniform_real_distribution<> dis(0.0, 1.0);
        for (size_t i = 0; i < individual.size(); ++i) {
            if (dis(gen) < mutation_rate) {
                std::uniform_int_distribution<> dis_idx(0, static_cast<int>(load_matrices[i].size()) - 1);
                individual[i] = dis_idx(gen);
            }
        }
    };

    // 初始化最佳解
    std::vector<int> best_solution;
    double best_score = std::numeric_limits<double>::infinity();
    double best_MUX_num = 0.0;

    int stagnation_count = 0;  // 追踪 best_score 连续相同的次数
    int max_stagnation = config.MAX_STAGNATION;  // 最大允许的连续相同次数
    double previous_best_score = 0.0;  // 用于记录前一代的 best_score

    // 主循环
    for (int generation = 0; generation < generations; ++generation) {
        // 评估种群
        std::vector<std::tuple<double, double, double>> scores(population_size);
        for (int i = 0; i < population_size; ++i) {
            scores[i] = evaluate(population[i]);
        }

        // 更新自适应交叉率和变异率
        double dynamic_crossover_rate = std::max(0.8, crossover_rate - 0.01 * (generation / static_cast<double>(generations)));
        double dynamic_mutation_rate = std::min(0.3, mutation_rate + 0.01 * (generation / static_cast<double>(generations)));

        // 提取成本、MUX_num和适应度
        std::vector<double> costs(population_size);
        std::vector<double> MUX_nums(population_size);
        std::vector<double> fitnesses(population_size);

        for (int i = 0; i < population_size; ++i) {
            costs[i] = std::get<0>(scores[i]);
            MUX_nums[i] = std::get<1>(scores[i]);
            fitnesses[i] = std::get<2>(scores[i]);
        }

        // 选择当前代中的最佳解（根据加权适应度值）
        auto min_it = std::min_element(fitnesses.begin(), fitnesses.end());
        int current_best_index = static_cast<int>(std::distance(fitnesses.begin(), min_it));
        auto current_best_solution = population[current_best_index];
        double current_best_score = costs[current_best_index];  // 仅保存cost
        double current_best_MUX_num = MUX_nums[current_best_index];

        if (fitnesses[current_best_index] < best_score) {
            best_solution = current_best_solution;
            best_score = current_best_score;
            best_MUX_num = current_best_MUX_num;
        }

        if (config.PRINT_VERILOG == true) {
            std::cout << "Generation " << generation << ", Best Cost: " << best_score << ", Best MUX: " << best_MUX_num << std::endl;
        }

        // 检查 best_score 是否连续相同
        if (best_score == previous_best_score) {
            stagnation_count += 1;
        } else {
            stagnation_count = 0;  // 重置计数器
            previous_best_score = best_score;  // 更新 previous_best_score
        }

        if (config.REPEAT_CHECK == true) {
            // 如果连续相同次数达到最大限制，则退出循环
            if (stagnation_count >= max_stagnation) {
                std::cout << "Best Cost稳定." << std::endl;
                break;
            }
        }

        // 按适应度排序种群
        std::vector<std::pair<double, std::vector<int>>> fitness_population;
        for (int i = 0; i < population_size; ++i) {
            fitness_population.emplace_back(fitnesses[i], population[i]);
        }
        std::sort(fitness_population.begin(), fitness_population.end(), [](const auto& a, const auto& b) {
            return a.first < b.first;
        });

        // 精英保留
        std::vector<std::vector<int>> new_population;
        for (int i = 0; i < elite_size; ++i) {
            new_population.push_back(fitness_population[i].second);
        }

        // 生成新种群
        while (static_cast<int>(new_population.size()) < population_size) {
            auto parents = select_parents(population, scores);
            auto child = crossover(parents[0], parents[1]);
            mutate(child);
            new_population.push_back(child);
        }

        population = new_population;
    }

    return std::make_tuple(best_solution, best_score, best_MUX_num);
}


std::vector<std::tuple<std::vector<int>, double, double>> run_ga_by_threads(const MatrixGroups& load_matrices) {
    // The number of threads
    ThreadPool threadPool(config.THREAD);

    auto GAConfig = dynamic_ga_params(load_matrices);

    // Results container for the best solutions
    std::vector<std::tuple<std::vector<int>, double, double>> results;

    // Mutex for results vector
    std::mutex results_mutex;

    // Enqueue tasks
    std::vector<std::future<void>> futures;

    for (int restart = 0; restart < config.NUM_RESTARTS; ++restart) {
        futures.emplace_back(std::async(std::launch::async, [restart, &results, &results_mutex, load_matrices, GAConfig]() {
            try {
                std::vector<int> best_solution;
                double best_score;
                double best_MUX_num;

                std::tie(best_solution, best_score, best_MUX_num) = genetic_algorithm(
                    load_matrices,
                    static_cast<int>(GAConfig.at("population_size")),
                    static_cast<int>(GAConfig.at("generations")),
                    GAConfig.at("crossover_rate"),
                    GAConfig.at("max_mutation_rate"),
                    GAConfig.at("elite_fraction")
                );
                {
                    std::lock_guard<std::mutex> lock(results_mutex);
                    results.emplace_back(best_solution, best_score, best_MUX_num);
                }

                std::cout << "Thread " << restart << " completed successfully with score: " << best_score << "\n";

            } catch (const std::exception& e) {
                std::cerr << "Thread " << restart << " encountered an error: " << e.what() << std::endl;
            }
        }));
    }

    for (auto& future : futures) {
        future.get(); // Wait for all threads to finish
    }

    return results;
}

MatrixGroup optimize_matrix_combination(const MatrixGroups& load_matrices) {
    // 获取所有矩阵的最大方阵尺寸
    int max_size = get_max_square_matrix_size(load_matrices);

    // 将每个 group 中的矩阵填充到最大方阵大小
    MatrixGroups padded_matrices;
    for (const auto& group : load_matrices) {
        MatrixGroup padded_group;
        for (const auto& matrix : group) {
            padded_group.push_back(pad_square_matrix(matrix, max_size));
        }
        padded_matrices.push_back(padded_group);
    }

    // 调用线程，启动GA算法
    std::vector<std::tuple<std::vector<int>, double, double>> results = run_ga_by_threads(load_matrices);
    // TODO 使用穷举+Jaccard 相似度获取 Top‑10 组合
    // std::vector<std::tuple<std::vector<int>, double, double>> results = exhaustive_top_k(padded_matrices, 10);

    // 按第二项（score）降序排序
    std::sort(results.begin(), results.end(),
              [](const auto& a, const auto& b) {
                  return std::get<1>(a) < std::get<1>(b);
              });

    // 保留前几个结果
    if (results.size() > config.RESULT_STORE_NUM) {
        results.resize(config.RESULT_STORE_NUM);
    }

    // 用于保存最优的 selected_matrices
    MatrixGroup bestSelectedMatrices;
    double bestMuxNum;

    // 假设已有的变量和函数声明
    // config、results、padded_matrices、bestSelectedMatrices、bestMuxNum 等

    if (config.STORE_ENABLE == false) {
        auto solution = std::get<0>(results[0]);
        bestMuxNum = std::get<2>(results[0]);
        for (size_t i = 0; i < solution.size(); ++i) {
            bestSelectedMatrices.push_back(padded_matrices[i][solution[i]]);
        }
    } else {
        // 确定存储文件夹名称
        std::string folderName = "data";
        // 判断文件夹是否存在，不存在则创建
        if (!std::filesystem::exists(folderName)) {
            std::filesystem::create_directories(folderName);
        }

        // 获取当前时间
        auto now = std::chrono::system_clock::now();
        std::time_t now_time = std::chrono::system_clock::to_time_t(now);
        // 格式化时间为字符串，格式为 YYYYMMDD_HHMMSS
        char timeStr[100];
        std::strftime(timeStr, sizeof(timeStr), "%Y%m%d_%H%M%S", std::localtime(&now_time));

        // 获取最佳得分，假设第一个结果为最佳解
        double best_score = std::get<1>(results[0]);

        // 生成新的文件名：当前时间 + "_" + best_score + ".txt"
        std::ostringstream filename;
        filename << folderName << "/" << timeStr << "_" << best_score << ".txt";

        // 打开文件
        std::ofstream outfile(filename.str());

        for (size_t rank = 0; rank < results.size(); ++rank) {
            const auto& result = results[rank];
            const auto& solution = std::get<0>(result);
            double score = std::get<1>(result);
            double muxNum = std::get<2>(result);

            // 生成 selected_matrices
            MatrixGroup selected_matrices;
            for (size_t i = 0; i < solution.size(); ++i) {
                selected_matrices.push_back(padded_matrices[i][solution[i]]);
            }

            // 合并生成 optimal_merged_matrix
            Matrix optimal_merged_matrix = selected_matrices[0];
            for (size_t i = 1; i < selected_matrices.size(); ++i) {
                optimal_merged_matrix = merge_matrices(optimal_merged_matrix, selected_matrices[i]);
            }

            Matrix A, B;
            std::tie(A, B) = merge_matrices_custom(selected_matrices, 0);

            // 写入文件
            outfile << "Solution " << rank+1 << ":\n";
            for (int idx : solution) outfile << idx << " ";
            outfile << "\nScore: " << score << "\nMuxNum: " << muxNum << "\n";

            outfile << "Selected Matrices:\n";
            for (const auto& matrix : selected_matrices) {
                print_matrix(matrix, outfile);
                outfile << "\n";
            }

            outfile << "Optimal Merged Matrix:\n";
            print_matrix(optimal_merged_matrix, outfile);
            outfile << "\n";

            outfile << "Matrix A:\n";
            print_matrix(A, outfile);
            outfile << "\n";

            outfile << "Matrix B:\n";
            print_matrix(B, outfile);
            outfile << "\n";

            outfile << "\n=============================\n";

            // 保存最优的 selected_matrices
            if (rank == 0) {
                bestSelectedMatrices = selected_matrices;
                bestMuxNum = muxNum;
            }
        }

        // 关闭文件
        outfile.close();
    }
    // std::ofstream file("mux_num.txt", std::ios::app);
    // if (file.is_open()) {
    //     file << bestMuxNum << "\n"; // 写入数据并换行
    // } else {
    //     std::cerr << "Error: Unable to open file " << "mux_num.txt" << std::endl;
    // }

    double bestScore = std::get<1>(results[0]);
    std::ofstream scoreFile("best_score.txt", std::ios::app);
    if (scoreFile.is_open()) {
        scoreFile << bestScore << "\n";  // 写入最佳score并换行
        scoreFile.close();
    } else {
        std::cerr << "Error: Unable to open file best_score.txt" << std::endl;
    }


    return bestSelectedMatrices;
}


} // namespace MatrixProject