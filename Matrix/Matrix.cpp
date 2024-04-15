//
// Created by konst on 06.09.2022.
//

#include <iostream>
#include <cmath>
#include <iomanip>
#include "Matrix.h"

Matrix::Matrix() : n(0), m(0), matrix(nullptr), k(1) { }

Matrix::Matrix(const Pair &dims)
: n(dims.first()), m(dims.second()), matrix(new double*[this->n]), k(1) {
    for (int i = 0; i < this->n; ++i)
        this->matrix[i] = new double[this->m] {0};
}

Matrix::Matrix(const Pair & dims, double scalar)
: n(dims.first()), m(dims.second()), matrix(new double*[this->n]), k(1) {
    for (int i = 0; i < this->n; ++i)
        this->matrix[i] = new double[this->m] {0};
    if (this->n != this->m) return;
    for (int i = 0; i < this->n; ++i)
        this->matrix[i][i] = scalar;
}

Matrix::Matrix(const Pair & dims, const double * arr)
: n(dims.first()), m(dims.second()), matrix(new double*[this->n]), k(1) {
    for (int i = 0; i < this->n; ++i)
        this->matrix[i] = new double[this->m] {0};
    if (this->n != this->m) return;
    for (int i = 0; i < this->n; ++i)
        this->matrix[i][i] = arr[i];
}

Matrix::Matrix(Matrix const &that)
    : n(that.n)
    , m(that.m)
    , matrix(new double*[that.n])
    , k(that.k) {
    for (int i = 0; i < n; ++i) {
        matrix[i] = new double[m];
        for (int j = 0; j < m; ++j) matrix[i][j] = that.matrix[i][j];
    }
}

Matrix::~Matrix() {
    for (int i = 0; i < this->n; ++i) {
        delete[] this->matrix[i];
    }
    delete[] this->matrix;
}

double Matrix::determinant(Matrix &M){
    if (M.n != M.m) throw "Matrix is not quad";
    if (M.n == 1) return M.matrix[0][0];
    double res = 0;

    return res;
}

Matrix Matrix::t() {
    Pair dim(this->m, this->n);
    Matrix M(dim);
    for(int i = 0; i < this->n; ++i)
        for(int j = 0; j < this->m; ++j)
            M.matrix[j][i] = this->matrix[i][j];
    return M;
}

Matrix Matrix::operator !() {
    return this->t();
}

Matrix Matrix::operator + (Matrix const &that) const {
    if (this->n == 0 || that.n == 0 || this->m == 0 || that.m == 0) throw std::string("Invalid matrix");
    if (this->n != that.n || this->m != that.m) throw std::string ("Matrices have different dimensions");
    Matrix M(Pair(this->n, this->m));
    for(int i = 0; i < this->n; ++i)
        for(int j = 0; j < this->m; ++j)
            M.matrix[i][j] = this->matrix[i][j] + that.matrix[i][j];
    return M;
}

Matrix Matrix::operator +() const {
    return *this;
}

void Matrix::operator +=(Matrix const &that) {
    if (this->n == 0 || that.n == 0 || this->m == 0 || that.m == 0) throw "Invalid matrix";
    if (this->n != that.n || this->m != that.m) throw "Matrices have different dimensions";
    for (int i = 0; i < this->n; ++i)
        for (int j = 0; j < this->m; ++j)
            this->matrix[i][j] += that.matrix[i][j];
}

Matrix Matrix::operator -(Matrix const &that) const {
    if (this->n == 0 || that.n == 0 || this->m == 0 || that.m == 0) throw "Invalid matrix";
    if (this->n != that.n || this->m != that.m) throw "Matrices have different dimensions";
    Matrix M(Pair(this->n, this->m));
    for(int i = 0; i < this->n; ++i)
        for(int j = 0; j < this->m; ++j)
            M.matrix[i][j] = this->matrix[i][j] - that.matrix[i][j];
    return M;
}

Matrix Matrix::operator -() const {
    Matrix M(Pair(this->n, this->m));
    for(int i = 0; i < this->n; ++i)
        for(int j = 0; j < this->m; ++j)
            M.matrix[i][j] = -this->matrix[i][j];
    return M;
}

void Matrix::operator -=(Matrix const &that) {
    if (this->n == 0 || that.n == 0 || this->m == 0 || that.m == 0) throw "Invalid matrix";
    if (this->n != that.n || this->m != that.m) throw "Matrices have different dimensions";
    for (int i = 0; i < this->n; ++i)
        for (int j = 0; j < this->m; ++j)
            this->matrix[i][j] -= that.matrix[i][j];
}

Matrix Matrix::operator *(Matrix const &that) const {
    if (this->m != that.n) throw;
    Matrix M(Pair(this->n, that.m));
    for (int i = 0; i < this->n; ++i)
        for (int j = 0; j < that.m; ++j)
            for (int k = 0; k < this->m; ++k)
                M.matrix[i][j] += this->matrix[i][k] * that.matrix[k][j];
    return M;
}

void Matrix::operator *=(Matrix const &that) {
    if (n == 0 || that.n == 0 || m == 0 || that.m == 0) throw std::string("Invalid matrix");
    if (m != that.n) throw std::string("Matrices have different dimensions");
    Matrix M(*this);
    for (int i = 0; i < this->n; ++i)
        for (int j = 0; j < that.m; ++j) {
            this->matrix[i][j] = 0;
            for (int k = 0; k < m; ++k)
                matrix[i][j] += M.matrix[i][k] * that.matrix[k][j];
        }
}

Matrix Matrix::operator /(Matrix const &that) const {
    // TODO: assertion
    return Matrix();
}

void Matrix::operator /=(Matrix const &that) {
    // TODO: assertion
}



Matrix Matrix::operator/(double scalar) const {
    Matrix N = *this;
    for (int i = 0; i < N.n; ++i)
        for (int j = 0; j < N.m; ++j)
            N.matrix[i][j] /= scalar;
    return N;
}

void Matrix::operator/=(double scalar) {
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            matrix[i][j] /= scalar;
}

Matrix & Matrix::operator =(Matrix const & that){
    if (n != that.n && m != that.m) throw;
    for (int i = 0; i < this->n; ++i)
        for (int j = 0; j < this->m; ++j)
            this->matrix[i][j] = that.matrix[i][j];
    return *this;
}

Row Matrix::operator[](unsigned index) {
    if (index < 0 || index >= this->n) throw "Invalid index";
    Row row(this->m, this->matrix[index]);
    return row;
}

Row Matrix::operator[](unsigned index) const {
    if (index < 0 || index >= this->n) throw "Invalid index";
    const Row row(this->m, this->matrix[index]);
    return row;
}

Column Matrix::operator()(unsigned index) {
    if (index < 0 || index >= this->m) throw "Invalid index";
    double ** arr = new double *[this->n];
    for (int i = 0; i < this->n; ++i)
        arr[i] = &this->matrix[i][index];
    Column column(this->n, arr);
    return column;
}

Matrix Matrix::operator()(unsigned index_i, unsigned index_j) {
    if(this->n == 0 || this->m == 0) throw "Invalid matrix";
    if (index_i < 0 || index_i >= this->n || index_j < 0 || index_j >= this->m) throw "Invalid index";
    --index_i; --index_j;
    Matrix M(Pair(this->n - 1, this->m - 1));
    for (int i = 0; i < M.n; ++i)
        for (int j = 0; j < M.m; ++j)
            M.matrix[i][j] = this->matrix[i + (i >= index_i)][j + (j >= index_j)];
    return M;
}

bool Matrix::operator ==(const Matrix &that) {
    if (this->n == 0 || that.n == 0) throw "Invalid matrix";
    if (this->n != that.n) throw "Matrices have different dimensions";
    bool res = true;
    for (int i = 0; res && i < this->n; ++i)
        for (int j = 0; res && j < this->m; ++j)
            res &= this->matrix[i][j] == that.matrix[i][j];
    return res;
}

bool Matrix::operator !=(const Matrix &that) {
    return !(*this == that);
}

std::ostream & operator <<(std::ostream & stream, const Matrix &M) {
    for (int i = 0; i < M.n; ++i) {
        for (int j = 0; j < M.m; ++j)
            stream << std::setw(5) << M.matrix[i][j];
        stream << std::endl;
    }
    return stream;
}

std::istream & operator >>(std::istream &stream, Matrix &M) {
    for (int i = 0; i < M.n; ++i)
        for (int j = 0; j < M.m; ++j)
            stream >> M.matrix[i][j];
    return stream;
}

Matrix operator*(double scalar, Matrix const &M) {
    Matrix N(M);
    for (int i = 0; i < N.n; ++i)
        for (int j = 0; j < N.m; ++j)
            N.matrix[i][j] *= scalar;
    return N;
}

double Matrix::determinantTriangle() {
    double res = this->k;
    for (int i = 0; i < this->n; ++i) {
        res *= this->matrix[i][i];
    }
    return res;
}



Matrix Matrix::gauss(Matrix & that) {
    if (that.n != this->n) throw "Matrices have different dimensions";
    
    Matrix res(*this);
    for (int i = 0; i < res.n; ++i) {
        if (res[i][i] == 0.) {
            int j = i + 1;
            while (j < res.m && res[j][i] == 0) j++;
            if (j == res.m) continue;
            that[i].swap(that[j]);
            res.k *= -1;
            res[i].swap(res[j]);
        }
        if (res[i][i] == 0.) continue;
        res.k *= res[i][i];
        that[i] /= res[i][i];
        res[i] /= res[i][i];
        for (int j = i + 1; j < res.m; ++j) {
            that[j] -= res[j][i] * that[i];
            res[j] -= res[j][i] * res[i];
        }
    }
    return res;
}

Matrix Matrix::gauss() {
    Matrix res(*this);
    for (int i = 0; i < res.n; ++i) {
        if (res[i][i] == 0.) {
            int j = i + 1;
            while (j < res.m && res[j][i] == 0) j++;
            if (j == res.m) continue;
            res.k *= -1;
            res[i].swap(res[j]);
        }
        if (res[i][i] == 0.) continue;
        res.k *= res[i][i];
        res[i] /= res[i][i];
        for (int j = i + 1; j < res.m; ++j) {
            res[j] -= res[j][i] * res[i];
        }
    }
    return res;
}

Matrix Matrix::reserveGauss(Matrix &that) {
    if (that.n != this->n) throw "Matrices have different dimensions";
    Matrix res(*this);

    for (int i = res.n - 1; i >= 0 ; --i)
        for (int j = i - 1; j >= 0; --j) {
            that[j] -= res[j][i] * that[i];
            res[j] -= res[j][i] * res[i];
        }

    return res;
}

Matrix Matrix::reserveGauss() {
    Matrix res(*this);

    for (int i = res.n - 1; i >= 0 ; --i)
        for (int j = i - 1; j >= 0; --j) {
            res[j] -= res[j][i] * res[i];
        }

    return res;
}

Matrix Matrix::createRotatingMatrix(int i, int j) {
    double a = this->matrix[i][i], b = this->matrix[j][i];
    double norm = std::sqrt(a * a + b * b);
    double cos = a / norm;
    double sin = b / norm;
    return createRotatingMatrix(i, j, cos, sin);
}

Matrix Matrix::createRotatingMatrix(int i, int j, double cos, double sin) {
    Matrix res(this->getDims(), 1);
    res.matrix[i][i] = res.matrix[j][j] = cos;
    res.matrix[j][i] = -(res.matrix[i][j] = sin);
    return res;
}

Pair Matrix::getDims() const {
    return Pair(n, m);
}

Matrix Matrix::rotating(Matrix &that) {
    if (that.n != this->n) throw "Matrices have different dimensions";
    Matrix res(*this);
    Pair dims = res.getDims();
    for (int i = 0; i < dims.first(); ++i)
        for (int j = i + 1; j < dims.first(); ++j) {
            that = res.createRotatingMatrix(i, j) * that;
            res = res.createRotatingMatrix(i, j) * res;
        }

    for (int i = 0; i < dims.first(); ++i) {
        res.k *= res[i][i];
        that[i] /= res[i][i];
        res[i] /= res[i][i];
    }
    return res;
}

Matrix Matrix::dopisat(Matrix & that) {
    Pair dims(this->n, this->m + that.m);
    Matrix res(dims);
    for (int i = 0; i < this->n; ++i) {
        for (int j = 0; j < this->m; ++j)
            res.matrix[i][j] = this->matrix[i][j];
        for (int j = this->m; j < dims.second(); ++j)
            res.matrix[i][j] = that.matrix[i][j - this->m];
    }
    return res;
}

Matrix Matrix::reflecting(Matrix &that) {
    Matrix res(*this);

    return res;
}

Matrix Matrix::createReflectingMatrix(int i) {
    Pair dimsW(this->n, 1);
    Matrix X(dimsW), Y(dimsW);
    double norm = 0;
    for (int j = i; j < dimsW.first(); ++j)
        X.matrix[j][0] = this->matrix[j][i];
    norm = sqrt((X.t() * X).matrix[0][0]);
    if (norm == 0) return Matrix(this->getDims(), 1);
    Y.matrix[i][0] = norm;
    Matrix W = X - Y;
    norm = (W.t() * W).matrix[0][0];
    if (norm == 0) return Matrix(this->getDims(), 1);
    return Matrix(this->getDims(), 1) - (2 / norm) * W * W.t();
}

Matrix Matrix::iteration(const Matrix &x_k, const Matrix &b, const Matrix &H) {
    Matrix r_k = *this * x_k - b;
    return x_k - H * r_k;
}

Matrix Matrix::iteration(const Matrix &x_k, const Matrix &b, double tau) {
    Matrix r_k = *this * x_k - b;
    return x_k - tau * r_k;
}

Matrix Matrix::GDA(const Matrix &b, double eps) {
    return GDA(b, b, eps);
}

Matrix Matrix::GDA(const Matrix &b, const Matrix &x_0, double eps) {
    Matrix x = x_0, r(b.getDims());
    int k = 0;
    double norm, normA;
    do {
        r = *this * x - b;
        norm = (r.t() * r).matrix[0][0];
        normA = ((*this * r).t() * r).matrix[0][0];
        if (normA == 0) return x;
        x = iteration(x, b, norm / normA);
        k++;
    } while (sqrt(norm) > eps);
    std::cout << k << std::endl;
    return x;
}

Matrix Matrix::GMRES(const Matrix &b, double eps) {
    return GMRES(b, b, eps);
}

Matrix Matrix::GMRES(const Matrix &b, const Matrix &x_0, double eps) {
    Matrix x = x_0, r(b.getDims()), thisR(r.getDims());
    int k = 0;
    double norm, normA, normAA;
    do {
        r = *this * x - b;
        thisR = *this * r;
        norm = (r.t() * r).matrix[0][0];
        normA = (thisR.t() * r).matrix[0][0];
        normAA = (thisR.t() * thisR).matrix[0][0];
        if (normAA == 0) return x;
        x = iteration(x, b, normA / normAA);
        k++;
    } while (sqrt(norm) > eps);
    std::cout << k << std::endl;
    return x;
}

Matrix Matrix::getVector(int i) {
    Matrix res(Pair(this->n, 1));
    for (int j = 0; j < this->n; ++j)
        res.matrix[j][0] = this->matrix[j][i];
    return res;
}

Matrix Matrix::polynomialIteration(const Matrix &x) const {
    Matrix res = x;
    double norm = sqrt((res.t() * res).matrix[0][0]);
    (res = *this * res) /= norm;
    return res;
}

std::pair<int, double> Matrix::eigenvectorWithMaxEigenvalue(double eps) {
    Matrix ort(Pair(this->n, 1));
    ort.matrix[0][0] = 1;
    return eigenvectorWithMaxEigenvalue(ort, eps);
}

std::pair<int, double> Matrix::eigenvectorWithMaxEigenvalue(const Matrix &x0, double eps) {
    Matrix x = x0, r(x0.getDims()), t(x0.getDims());
    double norm, normR; int k = 0;
    do {
        k++;
        x = polynomialIteration(x);
        norm = sqrt((x.t() * x).matrix[0][0]);
        r = *this * x - norm * x;
        normR = sqrt((r.t() * r).matrix[0][0]);
    } while (normR >= eps);
    return {k, norm};
}

std::pair<int, double> Matrix::eigenvectorWithMinEigenvalue(double eps) {
    Matrix ort(Pair(this->n, 1));
    ort.matrix[0][0] = 1;
    return eigenvectorWithMinEigenvalue(ort, eps);
}

double Matrix::norm1() {
    double res = 0, temp = 0;
    for (int i = 0; i < m; ++i) temp += matrix[0][i];
    for (int i = 1; i < n; ++i) {
        temp = 0;
        for (int j = 0; j < m; ++j) temp += matrix[i][j];
        res = std::max(res, temp);
    }
    return res;
}

std::pair<int, double> Matrix::eigenvectorWithMinEigenvalue(const Matrix &x0, double eps) {
    Matrix A = this->norm1() * Matrix(this->getDims(), 1) - *this;
    std::pair<int, double> pair = A.eigenvectorWithMaxEigenvalue(x0, eps);
    return {pair.first, this->norm1() - pair.second};
}

Matrix Matrix::Jacobi(double eps) {
    Matrix res = *this;

    int k = 0; double f;
    do {
        k++;
        res = res.JacobiRotate();
        f = 0;
        for (int a = 0; a < n; ++a)
            for (int b = 0; b < m; ++b)
                if (a == b) continue;
                else f += res.matrix[a][b] * res.matrix[a][b];
        f = sqrt(f);
    } while (f > eps);
    int numbers = -(int)std::log10(eps);
    std::cout << k << '\n' << std::fixed << std::setprecision(numbers) << res;
    return res;
}

Matrix Matrix::JacobiRotate() {
    Matrix res = *this;
    int p = 0, q = 1;
    double max = fabs(res.matrix[0][1]);
    for (int i = 0; i < n; ++i)
        for (int j = 1 + i; j < m; ++j)
            if (fabs(res.matrix[i][j]) > max) {
                max = fabs(res.matrix[i][j]);
                p = i; q = j;
                break;
            }
    if (p == 0 && q == 0) return *this;
    double aii = matrix[p][p], aij = matrix[p][q], ajj = matrix[q][q];
    double tau = (aii - ajj) / ( 2 * aij);
    double t = - tau + (tau > 0? 1: -1) * sqrt(tau * tau + 1);
    if (t < 1.e-8) t = 1;
    double cos = 1. / sqrt(t * t + 1);
    double sin = t * cos;
    Matrix Q = Matrix::createRotatingMatrix(p, q, cos, -sin);
    res = Q.t() * res * Q;
    return res;
}

double Matrix::dot(const Matrix &other) const {
    if (other.n != this->n or other.m != this->m)
        throw "Invalid matrix dimensions";
    double s = 0;
    for (int i = 0; i < other.n; i++)
        for (int j = 0; j < other.m; j++)
            s += other[i][j] * (*this)[i][j];
    return s;
}

double Matrix::norm2() const {
    return std::sqrt(this->dot(*this));
}

//void Matrix::Jacobi() {
//    Matrix p(getDims());
//    for (int i = 0; i < n; i++) p.matrix[i][i] = 1;
//    int k = 50 * (n * m);
//    for (int i = 0; i < k; i++) {
//        Rotation(p);
//        //system("pause");
//    }
//    std::cout << p;
//}
//
//void Matrix::maxInd(int &ind_i, int &ind_j) {
//    double max = matrix[0][1];
//    ind_i = 0;
//    ind_j = 1;
//    for(int i = 0; i < n - 1; i++)
//        for (int j = i + 1; j < m; j++) {
//            if (abs(matrix[i][j]) > max) {
//                max = abs(matrix[i][j]);
//                ind_i = i;
//                ind_j = j;
//            }
//        }
//}
//
//void Matrix::Rotation(Matrix& p) {
//    int a, b;
//    maxInd(a, b);
//
//    Matrix rot1 = *this;
//    (*this) = (rot1 * (*this)) * rot1.t();
//    p = rot1 * p;
//}
//
