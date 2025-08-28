#!/bin/bash
## ./local.sh
#    # 定义包含多个整数的数组
#    # Sa
#    # values=("362" "392" "473")
#    # Sb
#    #values=("39" "150" "192")
#    # Sc
##       values=("39" "45" "41" "47")
#    # Se
#    #    values=("162" "50" "26" "15")
##      values=("18" "27" "69" "111")
#    #   values=("64" "5:88" "19:109" "40:123")
#      #  values=("64" "88" "109" "123")
#    #values=("64")
##       values=("5" "19" "40")
#       values=("117" "118")
##      values=("18 : 27 : 69 : 111")
##      values=("2621" "7036" "17695" "28353") # Q15
#
#    # Sd
##    values=("8027" "7416:3072" "5676")
#      values=("2621" "7036" "17695" "28353") # Q15
#
#values=("22664" "86735" "181137" "289546" "393216" "474223" "518560")

#
##     初始化 values_1 和 values_2 数组
#     values_1=()
#     values_2=()
#
#     # 计算 log2
#     log2() {
#         echo "scale=10; l($1)/l(2)" | bc -l
#     }
#
#     # 计算所需位宽
#     get_bit_width() {
#         local num=$1
#         # 计算绝对值
#         abs_num=$(echo "$num" | awk '{if ($1<0) print -$1; else print $1}')
#
#         # 计算 log2(绝对值)
#         log2_value=$(log2 "$abs_num")
#
#         # 向上取整并加1
#         bit_width=$(echo "($log2_value + 0.999999999) / 1 + 1" | bc)
#
#         echo "$bit_width"
#     }
#
#     # 遍历 values 数组，处理每个元素
#     for value in "${values[@]}"; do
#         # 拆分数字对
#         IFS=":" read -r num1 num2 <<< "$value"
#
#         # 将数字分别添加到 values_1 和 values_2
#         values_1+=("$num1")
#         values_2+=("$num2")
#     done
#
#     # 输出 values_1 和 values_2
#     echo "values_1: ${values_1[@]}"
#     echo "values_2: ${values_2[@]}"
#
#     # 计算 values_1 和 values_2 的最大位宽
#     max_bit_width_1=0
#     max_bit_width_2=0
#
#     # 计算 values_1 中的最大位宽
#     for value in "${values_1[@]}"; do
#         bit_width=$(get_bit_width "$value")
#         if [ "$bit_width" -gt "$max_bit_width_1" ]; then
#             max_bit_width_1=$bit_width
#         fi
#     done
#
#     # 计算 values_2 中的最大位宽
#     for value in "${values_2[@]}"; do
#         bit_width=$(get_bit_width "$value")
#         if [ "$bit_width" -gt "$max_bit_width_2" ]; then
#             max_bit_width_2=$bit_width
#         fi
#     done
#
#     # 输出最大位宽
#     echo "values_1 最大位宽: $max_bit_width_1 位"
#     echo "values_2 最大位宽: $max_bit_width_2 位"
#
#    commands=("cadical" "z3")
#
#    # # 清除txt文件中的所有内容
#    truncate -s 0 matrix_output.txt
#    truncate -s 0 matrix_length.txt
#    # 找到系数集中加法器最大的数量
#    solver_name=cadical
#    timeout=300
#    threads=16
#    quiet=1
#    minimize_full_adders=0
#    allow_post_adder_right_shifts=1
#    allow_negative_coefficients=1
#    write_cnf_files=0
#    allow_coefficient_sign_inversion=0
#    min_num_adders=0
#    enumerate_all=0
#
#    # 清除txt文件中的所有内容
#    truncate -s 0 matrix_output.txt
#
#    #遍历数组中的每个元素
#    for i in "${!values[@]}"; do
#        value="${values[$i]}"
#        # 遍历每个求解器
#        for command in "${commands[@]}"; do
#            echo "$value, $command, $min_num_adders"
##            ./satmcm "$value" "$command" "$timeout" "$threads" "$quiet" "$minimize_full_adders"\
#            ./satmcm_sh "$value" "$command" "$timeout" "$threads" "$quiet" "$minimize_full_adders"\
#            "$allow_post_adder_right_shifts" "$allow_negative_coefficients" "$write_cnf_files"\
#            "$allow_coefficient_sign_inversion" "$min_num_adders" "1"
#        done
#        echo "done" >> matrix_output.txt
#    done
#    echo "max_bit_width: $max_bit_width"
#     values=("117")

     # stage 5
#     values=("19" "38" "57")

     # stage 6
#      values=("3")

#    values=("30709:10694")

#    values=("13903")
#    values=("1835")
#    0.4242861112477116
# W8
#    values=("6" "21" "44")
#    values=("71" "96" "116" "127")
#    values=("71" "6:96" "21:116" "44:127")
# W12
#    values=("89" "339" "708")
#    values=("1131" "1536" "1852" "2026")
#    values=("1131" "89:1536" "339:1852" "708:2026")
# W16
#  values=("1416" "5421" "11321")
#  values=("18097" "24576" "29639" "32410")
#  values=("18097" "1416:24576" "5421:29639" "11321:32410")
# W20
#    values=("22664" "86735" "181137")
#    values=("289546" "393216" "474223" "518560")
    values=("289546" "22664:393216" "86735:474223" "181137:518560")


    commands=("cadical" "z3")

   # # 清除txt文件中的所有内容
    truncate -s 0 matrix_output.txt
    truncate -s 0 matrix_length.txt
    # 找到系数集中加法器最大的数量
    solver_name=cadical
    timeout=30000
    threads=16
    quiet=1
    minimize_full_adders=0
    allow_post_adder_right_shifts=1
    allow_negative_coefficients=1
    write_cnf_files=0
    allow_coefficient_sign_inversion=0
    min_num_adders=1
    enumerate_all=0

    # 清除txt文件中的所有内容
    truncate -s 0 matrix_output.txt

    #遍历数组中的每个元素
    for i in "${!values[@]}"; do
        value="${values[$i]}"
        # 遍历每个求解器
        for command in "${commands[@]}"; do
            echo "$value, $command, $min_num_adders"
            ./satmcm "$value" "$command" "$timeout" "$threads" "$quiet" "$minimize_full_adders"\
            "$allow_post_adder_right_shifts" "$allow_negative_coefficients" "$write_cnf_files"\
            "$allow_coefficient_sign_inversion" "$min_num_adders" "0"
        done
        echo "done" >> matrix_output.txt
    done



    ../cmake-build-debug/MatrixProject -f "matrix_output.txt" --restarts 20 --repeat-check 1 --stagnation 30 --thread 8 --store-enable 1 --result-num 10 --fraction-enable 0 --transform 1  --bitwidth 20
