//
// Created by konst on 26.09.2022.
//

#ifndef TASK_ROW_H
#define TASK_ROW_H


class Row {
private:
    const unsigned m;
    double *row;
public:
    Row(unsigned m, double *row);
    Row(Row const & row);

    void operator +=(Row const &row);
    void operator -=(Row const &row);

    Row operator -();
    Row operator +(Row const & row);
    Row operator -(Row const & row);

    void operator = (Row const & row);

    Row operator *(double scalar);
    void operator *=(double scalar);
    Row operator /(double scalar);
    void operator /=(double scalar);

    void swap(Row const & row);

    friend Row operator *(double scalar, Row const & row);
    double& operator[](unsigned index);
    double operator[](unsigned index) const;
};


#endif //TASK_ROW_H
