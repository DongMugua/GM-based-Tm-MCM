# GM-Based Tm-MCM

## Project Description

TmCM-AM is a project for matrix analysis and optimization. It supports matrix reading, transformation, optimization, and generates corresponding Verilog code.

## Compile and run

### Dependencies

- CMake 3.8 or higher
- Compiler that supports C++17

### Compilation steps

```bash
mkdir build
cd build
cmake ..
make
```

### Command line explanation
The program supports the following command-line options, allowing for customized behavior:

- -t, --transform:
  - Enables or disables the matrix transformation feature (MATRIX_TRANSFORM).
  - Default is 1.
- -b, --bitwidth:
  - input width.
  - Default is 8.
- -c, --repeat-check:
  - Enables or disables the repeat-check feature (REPEAT_CHECK)
  - Default is 1.
- -s, --stagnation:
  - Sets the maximum stagnation limit (MAX_STAGNATION) to control convergence in the optimization algorithm, if REPEAT_CHECK is enable.
  - Default is 40.
- -v, --verilog:
  - Enables or disables the Verilog output generation (PRINT_VERILOG).
  - Default is 1 (enabled).
- -n, --restarts:
  - Specifies the number of restarts for the algorithm (NUM_RESTARTS).
  - Default is 1.
- -f, --file:
  - Specifies the path to the input matrix file (MATRIX_INPUT)
  - Default is "matrix_output.txt".
- -p, --thread:
  - Number of CPU threads.
- -q, --store-enable:
  - Store results into text enable (STORE_ENABLE).
- -r, --result-num:
  - Number of Results Storage in text, if STORE_ENABLE is enable.
- -a --fraction-enable:
  - 小数乘法器使能信号，1为开启，0为关闭。开启时需要将--transform设置为0。
  - Fractional multiplier enable signal, 1 for on, 0 for off. When on, you need to set --transform to 0.
- i --truncate-factor:
  - 在小数常数乘法器中，此参数的含义是所保留的最大移位位宽。……
  - 
- w --fraction-width:
  - 小数点的位置。在小数乘法中，一般将小数点
  - The position of the decimal point. In fractional multiplication, the decimal point is generally
- l --only-for-print:
  - for just print the HDL codes according to the adjacent matrices form .txt. The .txt can be found in the script/data after any run. The exanple command can be found in script/5_only_for_print_HDL.sh.

An example:
```aiignore
    ../cmake-build-debug/MatrixProject -f "matrix_output.txt" -n 1 -c 1 -p 12
    ../cmake-build-debug/MatrixProject -f "matrix_output.txt" -n 300 -c 1 -s 100 -p 12
    ../cmake-build-debug/MatrixProject -f "matrix_output.txt" --restarts 400 --repeat-check 1 --stagnation 100 --thread 12 --result-num 10
    ../cmake-build-debug/MatrixProject -f "matrix_output.txt" --restarts 20 --repeat-check 1 --stagnation 30 --thread 8 --store-enable 1 --result-num 10 --fraction-enable 0 --transform 1
```
