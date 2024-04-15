#include <iostream>
#include <iomanip>
#include <vector>
#include "cmath"
#include "Matrix/Matrix.h"

double a = 0.8;
double b = 1.2;

double f(double x, double y) {
    return 0.4 * sin(y) * (3 * x * x - 4);
}

double phi(double x, double y) {
    return x * x * sin(y);
}

bool isInnerAreaPoint(double x, double y) {
    return y <= (x + 2./3) && y <= (5./3 - x) &&
           x <= 1 && x >= 0 &&
           y <= 1 && y >= 0;
}

bool isEdgeAreaPoint(double x, double y) {
    return y == (x + 2./3) || y == (5./3 - x) ||
           x == 1 || x == 0 ||
           y == 1 || y == 0;
}

Matrix getPattern(int n, double h) {
    Pair dims(n);
    Matrix pattern(dims);
    for (int i = 0; i < n; i++)
        for(int j = 0; j < n; j++) {
            double x = h * i;
            double y = h * j;
            if (!isInnerAreaPoint(x, y)) pattern[i][j] = -1;
            else if (isEdgeAreaPoint(x, y)) pattern[i][j] = 0;
            else pattern[i][j] = 1;
        }
    for (int i = 1; i < n - 1; i++)
        for(int j = 1; j < n - 1; j++)
            if (pattern[i][j] > 0 && (pattern[i - 1][j] < 0 || pattern[i + 1][j] < 0
                                  || pattern[i][j - 1] < 0 || pattern[i][j + 1] < 0)) {
                pattern[i][j] = 0;
            }
    return pattern;
}

void cut(Matrix &matrix, const Matrix &pattern,
         double h, double (&edge)(double, double)) {
    unsigned int n = pattern.getDims().first();
    double x, y;
    for (int i = 0; i < n; i++)
        for(int j = 0; j < n; j++) {
            x = h * i;
            y = h * j;
            if (pattern[i][j] == 0) matrix[i][j] = edge(x, y);
            else if (pattern[i][j] != 0) matrix[i][j] = 0;
        }
}


Matrix step(const Matrix &Up, const Matrix &pattern, double h) {
    unsigned int n = Up.getDims().first();
    Matrix res = Up;
    for (int i = 1; i < n - 1; i++)
        for(int j = 1; j < n - 1; j++)
            if (pattern[i][j] > 0)
                res[i][j] = (- a*(Up[i-1][j] - 2 * Up[i][j] + Up[i+1][j]) - b*(Up[i][j-1] - 2 * Up[i][j] + Up[i][j+1]))/h/h;
    return res;
}

Matrix function(const Matrix &pattern, double h, double (&func)(double, double), bool edges = false) {
    Matrix res(pattern.getDims());
    int n = res.getDims().first();
    double temp;
    for (int i = 0; i < n; i++)
        for(int j = 0; j < n; j++) {
            temp = pattern[i][j];
            res[i][j] = (temp > 0 || temp == 0 && edges)? func(h*i, h*j): 0;
        }
    return res;
}

Matrix simple_iter(const Matrix &U, const Matrix &F, const Matrix &pattern, double h) {
    Matrix res = U;
    int n = res.getDims().first();
    for (int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            if (pattern[i][j] > 0)
                res[i][j] = (h * h * F[i][j]
                        + a * (U[i-1][j] + U[i+1][j])
                        + b * (U[i][j-1] + U[i][j+1])) / (2*a + 2*b);
    return res;
}

Matrix upper_relax(const Matrix &U, const Matrix &F, const Matrix &pattern, double h, double omega) {
    Matrix res = U;
    int n = res.getDims().first();
    for (int i = 0; i < n; i++)
        for(int j = 0; j < n; j++)
            if (pattern[i][j] > 0)
                res[i][j] = U[i][j] + omega * (h * h * F[i][j]
                        + a * (res[i-1][j] - 2 * U[i][j] + U[i+1][j])
                        + b * (res[i][j-1] - 2 * U[i][j] + U[i][j+1])) / (2*a + 2*b);
    return res;
}

Matrix method2(double h, double eps) {
    int n = 1 / h + 1;
    Matrix pattern = getPattern(n, h);
    Matrix F = function(pattern, h, f);
    Matrix Up = pattern;
    cut(Up, pattern, h, phi);
    Matrix U = Up, r(pattern.getDims()), energR(pattern.getDims()), rPrev = r;
    Matrix PHI = function(pattern, h, phi);
    int k = 0;
    double norm, normA, normAA;
    do {
        rPrev = r;
        r = step(U, pattern, h) - F;
        energR = step(r, pattern, h);
        normA = energR.dot(r);
        normAA = energR.norm2();
        normAA *= normAA;
        if (normA == 0) break;
        U -= (normA / normAA) * r;
        k++;
    } while ((rPrev - r).norm2() * h > eps);
    Matrix PHI_M = function(pattern, h, phi, true);
    double diff = (U - PHI_M).norm2() / PHI_M.norm2();
//    std::cout << " " << std::setw(4) << k << " " << std::setw(11) << diff;
    return U;
}

Matrix method3(double h, double eps, double omega) {
    int n = 1. / h + 1;
    Matrix pattern = getPattern(n, h);
    Matrix F = function(pattern, h, f);
    Matrix Up = pattern;
    Matrix r(pattern.getDims()), rp(pattern.getDims());
    cut(Up, pattern, h, phi);
    Matrix U = Up;
    int k = 0;
    do {
        Up = U;
        U = upper_relax(Up, F, pattern, h, omega);
        rp = r;
        r = step(U, pattern, h) - F;
        k++;
    } while((r - rp).norm2() / r.norm2() > eps);
    Matrix PHI_M = function(pattern, h, phi, true);
    double diff = (Up - PHI_M).norm2() / PHI_M.norm2();
    std::cout << " " << std::setw(4) << k << " " << std::setw(11) << diff;
    return U;
}

Matrix method4(double h, double eps) {
    int n = 1. / h + 1;
    Matrix pattern = getPattern(n, h);
    Matrix F = function(pattern, h, f);
    Matrix Up = pattern;
    Matrix r(pattern.getDims()), rp(pattern.getDims());
    cut(Up, pattern, h, phi);
    Matrix U = Up;
    int k = 0;
    do {
        Up = U;
        U = simple_iter(Up, F, pattern, h);
        rp = r;
        r = step(U, pattern, h) - F;
        k++;
    } while((r - rp).norm2() / r.norm2() > eps);
    Matrix PHI_M = function(pattern, h, phi, true);
    double diff = (Up - PHI_M).norm2() / PHI_M.norm2();
    std::cout << " " << std::setw(4) << k << " " << std::setw(11) << diff;
    return U;
}

int main() {
    std::vector<double> eps_s = {1.e-6, 1.e-7, 1.e-8};
    std::vector<double> hs = {1. / 10, 1. / 20, 1. / 40};
    std::cout << "| h \\ eps |";
    for (double eps: eps_s) std::cout << std::setw(17) << eps << " |";
    for (double h: hs) {
        std::cout << std::endl;
        std::cout << "| " << std::setw(7) << h <<" |";
        for (double eps: eps_s) {
//            Matrix A = method2(h, eps);
//            Matrix B = method3(h, eps, 1.1);
//            Matrix B = method3(h, eps, 2/(1 + sin(M_PI * h)));
            Matrix B = method4(h, eps);
            std::cout << " |";
        }
    }
}
/*
 h \ eps | 1e-6           | 1e-7           | 1e-8
   1/10  | 779.293208 96  | 779.300188 121 | 779.300844 145
         | 20.7029123 45  | 20.6995593 59  | 20.6991334 76
   1/20  | 3178.96801 277 | 3179.08524 374 | 3179.09719 472
         | 20.9716906 112 | 20.9091655 165 | 20.902345 223
   1/40  | 12776.5664 615 | 12778.4305 992 | 12778.6211 1376
         | 22.3773738 216 | 21.4702746 404 | 21.3701137 611
*/

/*
 h \ eps | 1e-6             | 1e-7             | 1e-8
  1/10   | 269  2.17106e-05 | 315  2.2133e-05  | 361  2.21751e-05
  1/20   | 1031 5.03827e-06 | 1217 6.011e-06   | 1402 6.11027e-06
  1/40   | 3909 1.34117e-06 | 4655 1.37289e-06 | 5400 1.58235e-06
*/

/*
 h \ eps | 1e-6             | 1e-7             | 1e-8
  1/10   | 151  4.67385e-05 | 190  1.60942e-05 | 229  2.15455e-05
  1/20   | 538  0.000309251 | 705  2.61916e-05 | 2422 0.000131162
  1/40   | 1749 0.001324    | 2422 0.000131162 | 3096 1.17456e-05
*/

/*
    Opti:
    | h \ eps |            1e-06 |            1e-07 |            1e-08 |
    |     0.1 |   36 2.21795e-05 |   39 2.21806e-05 |   44 2.21806e-05 |
    |    0.05 |   76 6.23561e-06 |   85 6.23611e-06 |   93 6.23614e-06 |
    |   0.025 |  163 1.61708e-06 |  178 1.61716e-06 |  192 1.61716e-06 |

    Non-opti:
    | h \ eps |            1e-06 |            1e-07 |            1e-08 |
    |     0.1 |  105 2.20105e-05 |  121 2.21641e-05 |  137  2.2179e-05 |
    |    0.05 |  410 5.64199e-06 |  478 6.17515e-06 |  547 6.23015e-06 |
    |   0.025 | 1505 8.43777e-07 | 1780 1.43484e-06 | 2055 1.59869e-06 |
 */