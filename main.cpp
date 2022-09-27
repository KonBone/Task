#include <iostream>
#include "Matrix/Matrix.h"

int main() {
    Pair dims(3);
    Matrix M(dims), E(dims, 1);
    std::cin >> M;
    Matrix Temp(M);

    Temp.gauss(E);
    Temp.reserveGauss(E);

    std::cout << M << "----- \n" << E << "----- \n" << M * E;
    return 0;
}
