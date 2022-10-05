#include <iostream>
#include "Matrix/Matrix.h"

int main() {
    Pair dims(2);
    Pair dimsV(2, 1);
    Matrix M(dims), V(dimsV), E(dims, 1);
    std::cin >> M >> V;
    double det = M.gauss().determinantTriangle();
    Matrix Temp = M.gauss(E);
    Temp.reserveGauss(E);
    Temp = M.gauss(V);
    Temp.reserveGauss(V);
    std::cout << "Result:\n" << V << "Inverse:\n" << E << "Determinant:\n" << det << "\nA * A^-1:\n" << M * E;
    return 0;
}
