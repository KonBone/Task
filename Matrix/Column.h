//
// Created by konst on 23.09.2022.
//

#ifndef TASK_COLUMN_H
#define TASK_COLUMN_H


class Column {
private:
    unsigned n;
    double ** col;
public:
    Column(unsigned n, double ** arr);
    ~Column();
    double & operator [](unsigned i);
    unsigned size() const;
};


#endif //TASK_COLUMN_H
