#include <iostream>
#include "Matrix/Matrix.h"

int main() {
    Pair dims(3);
    Pair dimsV(3, 1);
    Matrix M(dims), N(dims, 1), E(dims, 1.), V(dimsV);
    std::cin >> M >> V;
    Matrix Temp(M);

    Temp.gauss(V);

    Temp.reserveGauss(V);

    std::cout << V;
    if (E == Temp)
        std::cout << N;
    else
        std::cout << "Matrix has not inverse one";
    return 0;
}
