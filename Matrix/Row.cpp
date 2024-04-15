//
// Created by konst on 26.09.2022.
//

#include "Row.h"

Row::Row(unsigned int m, double *row): m(m), row(row) { }

Row::Row(const Row &row): m(row.m), row(new double[row.m]) {
    for (int i = 0; i < this->m; ++i)
        this->row[i] = row.row[i];
}

void Row::operator+=(Row const &row) {
    if (row.m != this->m) throw "Rows have different dimensions";
    for (int i = 0; i < row.m; ++i)
        this->row[i] += row.row[i];
}

void Row::operator-=(Row const &row) {
    if (row.m != this->m) throw "Rows have different dimensions";
    for (int i = 0; i < row.m; ++i)
        this->row[i] -= row.row[i];
}

Row Row::operator-() {
    Row res(*this);
    for (int i = 0; i < this->m; ++i)
        res.row[i] = -res.row[i];
    return res;
}

Row Row::operator*(double scalar) {
    Row res(*this);
    for (int i = 0; i < this->m; ++i)
        res.row[i] *= scalar;
    return res;
}

double& Row::operator[](unsigned int index) {
    if (index < 0 || index >= this->m) throw "Invalid index";
    return this->row[index];
}

double Row::operator[](unsigned int index) const {
    if (index < 0 || index >= this->m) throw "Invalid index";
    return this->row[index];
}

Row Row::operator+(const Row &row) {
    if (this->m != row.m) throw "Rows have different dimensions";
    Row res(*this);
    for (int i = 0; i < this->m; ++i)
        res.row[i] += row.row[i];
    return res;
}

Row Row::operator-(const Row &row) {
    if (this->m != row.m) throw "Rows have different dimensions";
    Row res(*this);
    for (int i = 0; i < this->m; ++i)
        res.row[i] -= row.row[i];
    return res;
}

void Row::operator=(const Row &row) {
    if (this->m != row.m) throw "Rows have different dimensions";
    for (int i = 0; i < this->m; ++i)
        this->row[i] = row.row[i];
}

Row Row::operator/(double scalar) {
    Row res(*this);
    for (int i = 0; i < this->m; ++i)
        res.row[i] *= scalar;
    return res;
}

void Row::operator*=(double scalar) {
    for (int i = 0; i < this->m; ++i)
        this->row[i] *= scalar;
}

void Row::operator/=(double scalar) {
    for (int i = 0; i < this->m; ++i)
        this->row[i] /= scalar;
}

Row operator*(double scalar, const Row &row) {
    Row res(row);
    for (int i = 0; i < row.m; ++i)
        res.row[i] *= scalar;
    return res;
}

void Row::swap(const Row &row) {
    if (this->m != row.m) throw "Rows have different dimensions";
    double t;
    for (int i = 0; i < this->m; ++i) {
        t = this->row[i];
        this->row[i] = row.row[i];
        row.row[i] = t;
    }
}
