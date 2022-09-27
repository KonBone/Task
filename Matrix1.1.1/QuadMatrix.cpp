//
// Created by konst on 06.09.2022.
//

#include "QuadMatrix.h"

QuadMatrix::QuadMatrix() : n(0), matr(nullptr) { }

QuadMatrix::QuadMatrix(int n)
        : n(n)
        , matr(new int*[n]) {
    if (n <= 0) throw "Invalid arguments";
    for (int i = 0; i < n; ++i)
        this->matr[i] = new int[n] {0};
}

QuadMatrix::QuadMatrix(int n, int scalar)
        : n(n)
        , matr(new int*[n]) {
    if (n <= 0) throw "Invalid arguments";
    for (int i = 0; i < n; ++i) {
        this->matr[i] = new int[n] {0};
        this->matr[i][i] = scalar;
    }
}

QuadMatrix::QuadMatrix(int n, const int * arr)
        : n(n)
        , matr(new int*[n]) {
    if (n <= 0) throw "Invalid arguments";
    for (int i = 0; i < n; ++i) {
        this->matr[i] = new int[n]{0};
        this->matr[i][i] = arr[i];
    }
}

QuadMatrix::QuadMatrix(QuadMatrix const &that)
        : n(that.n)
        , matr(new int*[that.n]) {
    for (int i = 0; i < n; ++i) {
        this->matr[i] = new int[this->n];
        for (int j = 0; j < this->n; ++j) this->matr[i][j] = that.matr[i][j];
    }
}

QuadMatrix::~QuadMatrix() {
    for (int i = 0; i < this->n; ++i) {
        delete[] this->matr[i];
    }
    delete[] this->matr;
}

int QuadMatrix::determinant(QuadMatrix &M){
    if (M.n == 0) throw "Invalid matrix";
    int res = 0;

    return res;
}


int * QuadMatrix::at(int i, int j) {
    if (i >= this->n || j >= this->n) throw "Invalid arguments";
    return &this->matr[i][j];
}

QuadMatrix QuadMatrix::t() {
    QuadMatrix M(this->n);
    for(int i = 0; i < this->n; ++i)
        for(int j = 0; j < this->n; ++j)
            M.matr[j][i] = this->matr[i][j];
    return M;
}

QuadMatrix QuadMatrix::operator !() {
    return this->t();
}

QuadMatrix QuadMatrix::operator + (QuadMatrix const &that) const {
    if (this->n == 0 || that.n == 0) throw "Invalid matrix";
    if (this->n != that.n) throw "Matrices are different dimensions";
    QuadMatrix M(this->n);
    for(int i = 0; i < this->n; ++i)
        for(int j = 0; j < this->n; ++j)
            M.matr[i][j] = this->matr[i][j] + that.matr[i][j];
    return M;
}

QuadMatrix QuadMatrix::operator +() const {
    return *this;
}

void QuadMatrix::operator +=(QuadMatrix const &that) {
    if (this->n == 0 || that.n == 0) throw "Invalid matrix";
    if (this->n != that.n) throw "Matrices are different dimensions";
    for (int i = 0; i < this->n; ++i)
        for (int j = 0; j < this->n; ++j)
            this->matr[i][j] += that.matr[i][j];
}

QuadMatrix QuadMatrix::operator -(QuadMatrix const &that) const {
    if (this->n == 0 || that.n == 0) throw "Invalid matrix";
    if (this->n != that.n) throw "Matrices are different dimensions";
    QuadMatrix M(this->n);
    for(int i = 0; i < this->n; ++i)
        for(int j = 0; j < this->n; ++j)
            M.matr[i][j] = this->matr[i][j] - that.matr[i][j];
    return M;
}

QuadMatrix QuadMatrix::operator -() const {
    QuadMatrix M(this->n);
    for(int i = 0; i < this->n; ++i)
        for(int j = 0; j < this->n; ++j)
            M.matr[i][j] = -this->matr[i][j];
    return M;
}

void QuadMatrix::operator -=(QuadMatrix const &that) {
    if (this->n == 0 || that.n == 0) throw "Invalid matrix";
    if (this->n != that.n) throw "Matrices are different dimensions";
    for (int i = 0; i < this->n; ++i)
        for (int j = 0; j < this->n; ++j)
            this->matr[i][j] -= that.matr[i][j];
}

QuadMatrix QuadMatrix::operator *(QuadMatrix const &that) const {
    if (this->n == 0 || that.n == 0) throw "Invalid matrix";
    if (this->n != that.n) throw "Matrices are different dimensions";
    QuadMatrix M(this->n);
    for (int i = 0; i < this->n; ++i)
        for (int j = 0; j < that.n; ++j)
            for (int k = 0; k < this->n; ++k)
                M.matr[i][j] += this->matr[i][k] * that.matr[k][j];
    return M;
}

void QuadMatrix::operator *=(QuadMatrix const &that) {
    if (this->n == 0 || that.n == 0) throw "Invalid matrix";
    if (this->n != that.n) throw "Matrices are different dimensions";
    QuadMatrix M(*this);
    for (int i = 0; i < this->n; ++i)
        for (int j = 0; j < that.n; ++j)
            for (int k = 0; k < this->n; ++k)
                this->matr[i][j] += M.matr[i][k] * that.matr[k][j];
}

QuadMatrix QuadMatrix::operator =(QuadMatrix &that){
    return that;
}

int * QuadMatrix::operator[](int index) {
    if (index < 0 || index >= this->n) throw "Invalid index";
    return this->matr[index];
}

Column QuadMatrix::operator()(int index) {
    if (index < 0 || index >= this->n) throw "Invalid index";
    int ** arr = new int*[this->n];
    for (int i = 0; i < this->n; ++i)
        arr[i] = &this->matr[i][index];
    Column column(this->n, arr);
    return column;
}

QuadMatrix QuadMatrix::operator()(int index_i, int index_j) {
    if(this->n == 0) throw "Invalid matrix";
    if (index_i < 0 || index_i >= this->n || index_j < 0 || index_j >= this->n) throw "Invalid index";
    --index_i; --index_j;
    QuadMatrix M(this->n - 1);
    for (int i = 0; i < M.n; ++i)
        for (int j = 0; j < M.n; ++j)
            M.matr[i][j] = this->matr[i + (i >= index_i)][j + (j >= index_j)];
    return M;
}

bool QuadMatrix::operator ==(const QuadMatrix &that) {
    if (this->n == 0 || that.n == 0) throw "Invalid matrix";
    if (this->n != that.n) throw "Matrices are different dimensions";
    bool res = true;
    for (int i = 0; res && i < this->n; ++i)
        for (int j = 0; res && j < this->n; ++j)
            res &= this->matr[i][j] == that.matr[i][j];
    return res;
}

bool QuadMatrix::operator !=(const QuadMatrix &that) {
    return !(*this == that);
}

std::ostream & operator <<(std::ostream & stream, const QuadMatrix &M) {
    for (int i = 0; i < M.n; ++i) {
        for (int j = 0; j < M.n; ++j)
            stream << M.matr[i][j] << " ";
        stream << std::endl;
    }
    return stream;
}

std::istream & operator >>(std::istream &stream, QuadMatrix &M) {
    for (int i = 0; i < M.n; ++i)
        for (int j = 0; j < M.n; ++j)
            stream >> M.matr[i][j];
    return stream;
}

QuadMatrix operator*(int scalar, QuadMatrix const &M) {
    QuadMatrix N(M);
    for (int i = 0; i < N.n; ++i)
        for (int j = 0; j < N.n; ++j)
            N.matr[i][j] *= scalar;
    return N;
}
