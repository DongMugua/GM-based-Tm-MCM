#ifndef MATRIXTOFRACTION_H
#define MATRIXTOFRACTION_H

#include "MatrixLoader.h"
#include <vector>

namespace MatrixProject {

    /**
     * @brief 对输入的矩阵进行“fractionTransform”操作，从第二列开始，
     *        对每列先找到列中所有Element中的最大值(绝对值最大)，
     *        然后对该列和对应行执行等价变换(列除以该最大值, 行乘以该最大值)。
     *        最终对每列执行该操作后返回变换后的矩阵。
     *
     * @param inputMatrix 输入的Matrix
     * @return 变换后的Matrix
     */
    Matrix fractionTransform(const Matrix& inputMatrix);

} // namespace MatrixProject
#endif //MATRIXTOFRACTION_H
