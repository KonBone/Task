//
// Created by konst on 23.09.2022.
//

#ifndef TASK_COLUMN_H
#define TASK_COLUMN_H


class Column {
private:
    int n;
    int ** col;
public:
    Column(int n, int ** arr);
    ~Column();
    int & operator [](int i);
};


#endif //TASK_COLUMN_H
