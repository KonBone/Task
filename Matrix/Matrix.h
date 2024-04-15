//
// Created by konst on 06.09.2022.
//

#ifndef TASK_MATRIX_H
#define TASK_MATRIX_H
#include <string>
#include <istream>
#include "Row.h"
#include "Column.h"
#include "Pair.h"

class Matrix {
private:
    const unsigned n, m;
    double **matrix;
    double k;
public:
    Matrix();
    Matrix(const Pair &dims);
    Matrix(const Pair &dims, double scalar);
    Matrix(const Pair &dims, const double *arr);

    Matrix(Matrix const & that);
    ~Matrix();

    static double determinant(Matrix &M);
    double determinantTriangle();

    Matrix gauss(Matrix & that);
    Matrix gauss();

    Matrix rotating(Matrix & that);
    Matrix rotating();
    Matrix reflecting(Matrix & that);

    Matrix reserveGauss(Matrix & that);
    Matrix reserveGauss();

    Matrix t();
    Matrix operator !();

    Matrix operator + (Matrix const & that) const;
    Matrix operator + () const;
    void operator += (Matrix const & that);

    Matrix operator - (Matrix const & that) const;
    Matrix operator - () const;
    void operator -= (Matrix const & that);

    Matrix operator * (Matrix const & that) const;
    void operator *= (Matrix const & that);

    Matrix operator / (Matrix const & that) const;
    Matrix operator / (double scalar) const;
    void operator /= (Matrix const & that);
    void operator /= (double scalar);

    Matrix & operator = (Matrix const & that);

    Row operator [](unsigned index);
    Row operator [](unsigned index) const;
    Column operator ()(unsigned index);
    Matrix operator ()(unsigned i, unsigned j);

    bool operator ==(const Matrix & that);
    bool operator !=(const Matrix & that);

    friend std::ostream & operator <<(std::ostream & stream, const Matrix &M);
    friend std::istream & operator >>(std::istream & stream, Matrix &M);

    friend Matrix operator *(double scalar, Matrix const &M);

    Pair getDims() const;
    Matrix createRotatingMatrix(int i, int j);
    Matrix createRotatingMatrix(int i, int j, double cos, double sin);
    Matrix createReflectingMatrix(int i);
    Matrix dopisat(Matrix & that);

    Matrix iteration(const Matrix &x_k, const Matrix &b, const Matrix &H);
    Matrix iteration(const Matrix &x_k, const Matrix &b, double tau);

    Matrix GDA(const Matrix &b, double eps = 1.e-6);
    Matrix GDA(const Matrix &b, const Matrix &x_0, double eps = 1.e-6);

    Matrix GMRES(const Matrix &b, double eps = 1.e-6);
    Matrix GMRES(const Matrix &b, const Matrix &x_0, double eps = 1.e-6);

    double norm1();
    double norm2() const;
    double dot(const Matrix &other) const;
    Matrix getVector(int i);
    Matrix polynomialIteration(const Matrix &x) const;
    std::pair<int, double> eigenvectorWithMaxEigenvalue(double eps = 1.e-6);
    std::pair<int, double> eigenvectorWithMaxEigenvalue(const Matrix &x0, double eps = 1.e-6);
    std::pair<int, double> eigenvectorWithMinEigenvalue(double eps = 1.e-6);
    std::pair<int, double> eigenvectorWithMinEigenvalue(const Matrix &x0, double eps = 1.e-6);

    Matrix Jacobi(double eps = 1.e-6);
    Matrix JacobiRotate();
//    void Jacobi();
//    void maxInd(int &ind_i, int &ind_j);
//    void Rotation(Matrix& p);
};



#endif //TASK_MATRIX_H
