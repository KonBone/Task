//
// Created by konst on 06.09.2022.
//

#include <iostream>
#include <cmath>
#include "Matrix.h"

Matrix::Matrix() : n(0), m(0), matrix(nullptr), k(1) { }

Matrix::Matrix(Pair dims)
: n(dims.first()), m(dims.second()), matrix(new double*[this->n]), k(1) {
    for (int i = 0; i < this->n; ++i)
        this->matrix[i] = new double[this->m] {0};
}

Matrix::Matrix(Pair dims, double scalar)
: n(dims.first()), m(dims.second()), matrix(new double*[this->n]), k(1) {
    for (int i = 0; i < this->n; ++i)
        this->matrix[i] = new double[this->m] {0};
    if (this->n != this->m) return;
    for (int i = 0; i < this->n; ++i)
        this->matrix[i][i] = scalar;
}

Matrix::Matrix(Pair dims, const double * arr)
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
        this->matrix[i] = new double[this->m];
        for (int j = 0; j < this->m; ++j) this->matrix[i][j] = that.matrix[i][j];
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
    Matrix M(this->m, this->n);
    for(int i = 0; i < this->n; ++i)
        for(int j = 0; j < this->m; ++j)
            M.matrix[j][i] = this->matrix[i][j];
    return M;
}

Matrix Matrix::operator !() {
    return this->t();
}

Matrix Matrix::operator + (Matrix const &that) const {
    if (this->n == 0 || that.n == 0 || this->m == 0 || that.m == 0) throw "Invalid matrix";
    if (this->n != that.n || this->m != that.m) throw "Matrices have different dimensions";
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
    // TODO: assertion
    Matrix M(Pair(this->n, that.m));
    for (int i = 0; i < this->n; ++i)
        for (int j = 0; j < that.m; ++j)
            for (int k = 0; k < this->m; ++k)
                M.matrix[i][j] += this->matrix[i][k] * that.matrix[k][j];
    return M;
}

void Matrix::operator *=(Matrix const &that) {
    if (this->n == 0 || that.n == 0 || this->m == 0 || that.m == 0) throw "Invalid matrix";
    if (this->n != that.n) throw "Matrices have different dimensions";
    Matrix M(*this);
    for (int i = 0; i < this->n; ++i)
        for (int j = 0; j < that.m; ++j)
            for (int k = 0; k < this->m; ++k)
                this->matrix[i][j] += M.matrix[i][k] * that.matrix[k][j];
}

Matrix Matrix::operator /(Matrix const &that) const {
    // TODO: assertion
    return Matrix();
}

void Matrix::operator /=(Matrix const &that) {
    // TODO: assertion
}

Matrix & Matrix::operator =(Matrix const & that){
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
            stream << M.matrix[i][j] << " ";
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
    Matrix res(this->getDims(), 1);
    double a = this->matrix[i][i], b = this->matrix[j][i];
    double norm = std::sqrt(a * a + b * b);
    double cos = a / norm;
    double sin = -b / norm;
    res.matrix[i][i] = res.matrix[j][j] = cos;
    res.matrix[j][i] = sin;
    res.matrix[i][j] = -sin;
    return res;
}

Pair Matrix::getDims() const {
    return Pair(this->n);
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
