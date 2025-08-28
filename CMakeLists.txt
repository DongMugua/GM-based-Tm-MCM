# 指定最低 CMake 版本
cmake_minimum_required(VERSION 3.10)
project(MatrixProject)

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 如果未设置编译类型，默认为 Debug
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# 添加 include 目录
include_directories(include)

# 源文件列表
set(SOURCES
        src/Main.cpp
        src/MatrixLoader.cpp
        src/MatrixTools.cpp
        src/MatrixTools_for_HDL.cpp
        src/MatrixCost.cpp
        src/MatrixAB.cpp
        src/MatrixOptimizer.cpp
        src/MatrixTransform.cpp
        src/MatrixToHDL.cpp
        src/thread.cpp
        src/VerilogGenerator.cpp
        src/VerilogGenerator_fraction.cpp
        src/MatrixToFraction.cpp
)

# 配置编译器优化选项
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # GNU 和 Clang 的优化选项
    set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG -flto")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    # MSVC 的优化选项
    set(CMAKE_CXX_FLAGS_DEBUG "/Od /Zi")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /DNDEBUG")
endif()

#if(APPLE)
#    set(CMAKE_OSX_SYSROOT "/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk")
#endif()

# 启用链接时优化 (LTO)（可选）
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)

# 添加可执行文件
add_executable(MatrixProject ${SOURCES})

# 链接库（如果需要其他数学库可以在此添加）
# target_link_libraries(MatrixProject PRIVATE pthread)