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
    const unsigned n = 0, m = 0;
    double **matrix;
public:
    Matrix();
    Matrix(Pair dims);
    Matrix(Pair dims, double scalar);
    Matrix(Pair dims, const double *arr);

    Matrix(Matrix const & that);
    ~Matrix();

    static double determinant(Matrix &M);

    void gauss(Matrix & that);
    void reserveGauss(Matrix & that);

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
    void operator /= (Matrix const & that);

    Matrix & operator = (Matrix const & that);

    Row operator [](unsigned index);
    Column operator ()(unsigned index);
    Matrix operator ()(unsigned i, unsigned j);

    bool operator ==(const Matrix & that);
    bool operator !=(const Matrix & that);

    friend std::ostream & operator <<(std::ostream & stream, const Matrix &M);
    friend std::istream & operator >>(std::istream & stream, Matrix &M);

    friend Matrix operator *(double scalar, Matrix const &M);
};



#endif //TASK_MATRIX_H
