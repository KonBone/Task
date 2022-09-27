//
// Created by konst on 23.09.2022.
//

#include "Column.h"

Column::Column(int n, int **arr): n(n), col(arr) { }

Column::~Column() {
    delete[] this->col;
}

int &Column::operator[](int i) {
    return *this->col[i];
}