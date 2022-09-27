//
// Created by konst on 21.09.2022.
//

#ifndef TASK_QUADMATRIX_H
#define TASK_QUADMATRIX_H

#include <istream>
#include "Column.h"

class QuadMatrix {
private:
    const int n = 0;
    int **matr;
public:
    QuadMatrix();
    QuadMatrix(int n);
    QuadMatrix(int n, int scalar);
    QuadMatrix(int n, const int *arr);

    QuadMatrix(QuadMatrix const & that);
    ~QuadMatrix();

    static int determinant(QuadMatrix &M);

    int * at(int i, int j);

    QuadMatrix t();
    QuadMatrix operator !();

    QuadMatrix operator + (QuadMatrix const & that) const;
    QuadMatrix operator + () const;
    void operator += (QuadMatrix const & that);

    QuadMatrix operator - (QuadMatrix const & that) const;
    QuadMatrix operator - () const;
    void operator -= (QuadMatrix const & that);

    QuadMatrix operator * (QuadMatrix const & that) const;
    void operator *= (QuadMatrix const & that);

    QuadMatrix operator = (QuadMatrix & that);

    int * operator [](int index);
    Column operator ()(int index);
    QuadMatrix operator ()(int i, int j);

    bool operator ==(const QuadMatrix & that);
    bool operator !=(const QuadMatrix & that);

    friend std::ostream & operator <<(std::ostream & stream, const QuadMatrix &M);
    friend std::istream & operator >>(std::istream & stream, QuadMatrix &M);

    friend QuadMatrix operator *(int scalar, QuadMatrix const &M);
};


#endif //TASK_QUADMATRIX_H
