#ifndef CONFIG_H
#define CONFIG_H

#include "MatrixLoader.h"

namespace MatrixProject {
    struct Config {
        bool MATRIX_TRANSFORM = true;
        bool REPEAT_CHECK = true;
        int MAX_STAGNATION = 20;
        int BIT_WIDTH = 8;
        bool PRINT_VERILOG = true;
        int NUM_RESTARTS = 1;
        int THREAD = 4;
        bool STORE_ENABLE = true;
        bool FRACTION_INPUT = false;
        int FRACTION_WIDTH = 0;
        int RESULT_STORE_NUM = 1;
        int TRUNCATE_FACTOR = 16; //truncate-factor
        int OUT_WIDTH = BIT_WIDTH;
        std::string MATRIX_INPUT = "../script/matrix_output.txt";
        bool ONLY_FOR_PRINT = false;
    };

    // 将配置参数作为全局变量或在 main 函数中传递
    extern Config config;

}

#endif //CONFIG_H
