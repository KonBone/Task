//
// Created by konst on 27.09.2022.
//

#ifndef TASK_PAIR_H
#define TASK_PAIR_H


class Pair {
private:
    unsigned n, m;
public:
    Pair();
    Pair(unsigned n);
    Pair(unsigned n, unsigned m);

    unsigned first() const;
    unsigned second() const;

    unsigned operator[](unsigned i) const;
};


#endif //TASK_PAIR_H
