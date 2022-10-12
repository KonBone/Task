#include <iostream>
#include "Matrix/Matrix.h"

int main() {
    Pair dims(4), dimsV(4, 1);
    Matrix M(dims), V(dimsV), E(dims, 1);
    std::cin >> M >> V;
    Matrix right = V.dopisat(E);
    Matrix temp = M.rotating(right);
    temp = temp.reserveGauss(right);
    std::cout << right << temp.determinantTriangle();
    return 0;
}
/*
1 1 1 1
1 9.9 99.8 999.7
1 99.8 9999.9 999999.8
1 999.7 999999.8 999999999.9
10
4319
4030199.5
4003001999.4
 */
