#include "MatrixLoader.h"
#include "MatrixOptimizer.h"
#include "MatrixTools.h"
#include "MatrixAB.h"
#include "MatrixToHDL.h"
#include "MatrixTransform.h"
#include "VerilogGenerator.h"
#include "VerilogGenerator_fraction.h"
#include "MatrixToFraction.h"
#include "Config.h"
#include <iostream>
#include <chrono>
#include <unistd.h>
#include <cstdlib>
#include <algorithm>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <ctime>
#include <string>

namespace MatrixProject {
    Config config;  // 定义全局的 config 变量
}

void saveVerilogToFile(const std::string& verilog_code, const std::string& filename = "output_verilog.txt") {
    if (verilog_code.empty()) {
        std::cerr << "Error: Verilog code is empty, file not created!" << std::endl;
        return;
    }

    std::ofstream file(filename, std::ios::out);
    if (!file) {
        std::cerr << "Error: Unable to open file " << filename << " for writing!" << std::endl;
        return;
    }

    file << verilog_code;
    std::cout << "Verilog code saved to: " << filename << std::endl;
}

int main(int argc, char* argv[]) {
      using namespace MatrixProject;

    // 定义长选项
    static struct option long_options[] = {
        {"transform", required_argument, 0, 't'},
        {"bitwidth", required_argument, 0, 'b'},
        {"repeat-check", required_argument, 0, 'c'},
        {"stagnation", required_argument, 0, 's'},
        {"verilog", required_argument, 0, 'v'},
        {"restarts", required_argument, 0, 'n'},
        {"file", required_argument, 0, 'f'},
        {"thread", required_argument, 0, 'p'},
        {"result-num", required_argument, 0, 'r'},
        {"store-enable", required_argument, 0, 'q'},
        {"fraction-enable", required_argument, 0, 'a'},
        {"truncate-factor", required_argument, 0, 'i'},
        {"fraction-width", required_argument, 0, 'w'},
        {"out-width", required_argument, 0, 'o'},
        {"only-for-print", required_argument, 0, 'l'},
        {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;

    while ((opt = getopt_long(argc, argv, "t:b:c:s:v:n:f:p:r:q:", long_options, &option_index)) != -1) {
        switch (opt) {
            case 't':
                config.MATRIX_TRANSFORM = std::stoi(optarg) != 0;
                break;
            case 'b':
                config.BIT_WIDTH = std::stoi(optarg);
                break;
            case 'c':
                config.REPEAT_CHECK = std::stoi(optarg) != 0;
                break;
            case 's':
                config.MAX_STAGNATION = std::stoi(optarg);
                break;
            case 'v':
                config.PRINT_VERILOG = std::stoi(optarg) != 0;
                break;
            case 'n':
                config.NUM_RESTARTS = std::stoi(optarg);
                break;
            case 'f':
                config.MATRIX_INPUT = optarg;
                break;
            case 'p':
                config.THREAD = std::stoi(optarg);
                break;
            case 'r':
                config.RESULT_STORE_NUM = std::stoi(optarg);
                break;
            case 'q':
                config.STORE_ENABLE = std::stoi(optarg) != 0;
                break;
            case 'a':
                config.FRACTION_INPUT = std::stoi(optarg);
                break;
            case 'w':
                config.FRACTION_WIDTH = std::stoi(optarg);
                break;
            case 'i':
                config.TRUNCATE_FACTOR = std::stoi(optarg);
                break;
            case 'o':
                config.OUT_WIDTH = std::stoi(optarg);
                break;
            case 'l':
                config.ONLY_FOR_PRINT = std::stoi(optarg);
                break;
            default:
                std::cerr << "Usage: " << argv[0] << " [options]\n";
                std::cerr << "Options:\n";
                std::cerr << "  -t, --transform     MATRIX_TRANSFORM (default 1)\n";
                std::cerr << "  -b, --bitwidth      BIT_WIDTH (default 8)\n";
                std::cerr << "  -c, --repeat-check  REPEAT_CHECK (default 1)\n";
                std::cerr << "  -s, --stagnation    MAX_STAGNATION (default 40)\n";
                std::cerr << "  -v, --verilog       PRINT_VERILOG (default 1)\n";
                std::cerr << "  -n, --restarts      NUM_RESTARTS (default 20)\n";
                std::cerr << "  -f, --file          MATRIX_INPUT file (default \"matrix_output.txt\")\n";
                std::cerr << "  -p, --thread NUM    THREAD count (default 1)\n";
                std::cerr << "  -r, --result-store  NUM Result Store num count (default 1)\n";
                std::cerr << "  -q, --store-enable  STORE_ENABLE (default 1)\n";
                std::cerr << "  -a, --fraction-enable   fraction input enable (default 0)\n";
                std::cerr << "  -i, --truncate-factor  TRUNCATE_FACTOR (default 0)\n";
                std::cerr << "  -w, --fraction-width  fraction-width (default 0)\n";
                std::cerr << "  -o, --out-width  ouput-width (default 8)\n";
                exit(1);
        }
    }

    if (!config.ONLY_FOR_PRINT) {

      // ***************************读取矩阵******************************************
        MatrixGroups load_matrices;

        load_matrices = read_matrix_from_txt(config.MATRIX_INPUT);

            // ***************************矩阵移除一些数据***************
        // if (config.DEBUG_REMOVE) {
        //     // 获取第二组数据的引用
        //     MatrixGroup& second_group = load_matrices[1];
        //     // 筛选出符合条件的矩阵
        //     MatrixGroup new_second_group;
        //     for (const auto& matrix : second_group) {
        //         if (/* 判断条件 */1) {
        //             new_second_group.push_back(matrix);
        //         }
        //     }
        //     second_group = new_second_group;
        // }

        // ***************************变换矩阵******************************************
        // 矩阵去重
        load_matrices = remove_duplicate_matrices(load_matrices);

        if (config.FRACTION_INPUT) {
            for (auto &matrix_group : load_matrices) {
                for (auto &mat : matrix_group) {
                    print_matrix(mat);
                    mat = fractionTransform(mat);
                    print_matrix(mat);
                    std::cout << std::endl;
                }
            }
        }

        std::vector<MatrixGroup> augmented_matrix_groups;
        if (config.MATRIX_TRANSFORM && !config.FRACTION_INPUT) {
            for (const auto& matrix_group : load_matrices) {
                augmented_matrix_groups.push_back(matrix_transform(matrix_group, 0.125));
            }
        } else {
            augmented_matrix_groups = load_matrices;
        }

        // Analyze the maximum matrix size
        int max_size = 0;
        for (const auto& group : augmented_matrix_groups) {
            for (const auto& matrix : group) {
                if (matrix.size() == matrix[0].size()) {
                    max_size = std::max(max_size, static_cast<int>(matrix.size()));
                }
            }
        }

        // 矩阵填充
        std::vector<MatrixGroup> padded_matrix = pad_all_groups(augmented_matrix_groups);

        // 确保所有元素都是列表格式
        for (auto& group : padded_matrix) {
            for (auto& matrix : group) {
                matrix = ensure_list_format(matrix);
            }
        }
        print_matrix(padded_matrix);
        // 矩阵去重
        augmented_matrix_groups = remove_duplicate_matrices(padded_matrix);

        print_matrix(augmented_matrix_groups);

        // Matrix matrix392 = {
        //     {std::vector<double>{0.0}, std::vector<double>{8.0, -32.0}, std::vector<double>{0.0}, std::vector<double>{8.0}},
        //     {std::vector<double>{0.0}, std::vector<double>{0.0}, std::vector<double>{0.0}, std::vector<double>{-16.0}},
        //     {std::vector<double>{0.0}, std::vector<double>{0.0}, std::vector<double>{0.0}, std::vector<double>{0.0}},
        //     {std::vector<double>{0.0}, std::vector<double>{0.0}, std::vector<double>{0.0}, std::vector<double>{1.0}}
        // };
        //
        // // 搜索 matrix392
        // // 测试是否存在
        // if (matrix_exists(matrix392, augmented_matrix_groups[1])) {
        //     std::cout << "Matrix exists in the group!" << std::endl;
        // } else {
        //     std::cout << "Matrix does not exist in the group." << std::endl;
        // }


      // ***************************矩阵从 group 中选择*****************************************
        auto start_time = std::chrono::high_resolution_clock::now();

        print_matrix(augmented_matrix_groups);

        const MatrixGroup selected_matrices = optimize_matrix_combination(augmented_matrix_groups);

        auto end_time = std::chrono::high_resolution_clock::now();
        auto execution_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();

        std::cout << "Select算法执行时间: " << execution_time << " s" << std::endl;

        // ***************************合并矩阵******************************************
        // 实现 merge_matrices_custom 函数
        Matrix A, B;
        std::tie(A, B) = merge_matrices_custom(selected_matrices, 1);

        // ***************************生成Verilog代码******************************************
        if (config.PRINT_VERILOG) {
            // 生成 Verilog 代码
            if (config.FRACTION_INPUT) {
            std::string verilog_code = generate_verilog_with_matrices_fraction(selected_matrices, A, B, config.BIT_WIDTH);
            std::cout << verilog_code << std::endl;
            saveVerilogToFile(verilog_code);
            }else{
                std::string verilog_code = generate_verilog_with_matrices(selected_matrices, A, B, config.BIT_WIDTH, config.OUT_WIDTH);
            std::cout << verilog_code << std::endl;
            saveVerilogToFile(verilog_code);
            }

        }

        return 0;

    } else {
        std::vector<Matrix> selected_matrices;
        Matrix A, B;

        // 调用我们编写的解析函数
        if (!parse_matrices_from_file(config.MATRIX_INPUT, selected_matrices, A, B)) {
            std::cerr << "解析矩阵文件失败！" << std::endl;
            return -1;
        }

        if (config.FRACTION_INPUT) {
            std::string verilog_code = generate_verilog_with_matrices_fraction(selected_matrices, A, B, config.BIT_WIDTH);
            std::cout << verilog_code << std::endl;
            saveVerilogToFile(verilog_code);
        }else{
            std::string verilog_code = generate_verilog_with_matrices(selected_matrices, A, B, config.BIT_WIDTH, config.OUT_WIDTH);
            std::cout << verilog_code << std::endl;
            saveVerilogToFile(verilog_code);
        }

        return 0;
    }
}