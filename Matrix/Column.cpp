//
// Created by konst on 23.09.2022.
//

#include "Column.h"

Column::Column(unsigned n, double **arr): n(n), col(arr) { }

Column::~Column() {
    delete[] this->col;
}

double &Column::operator[](unsigned i) {
    return *this->col[i];
}

unsigned Column::size() const {
    return n;
}
