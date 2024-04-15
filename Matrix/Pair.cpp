//
// Created by konst on 27.09.2022.
//

#include "Pair.h"

Pair::Pair() : n(0), m(0) { }

Pair::Pair(unsigned int n) : n(n), m(n) { }

Pair::Pair(unsigned int n, unsigned int m) : n(n), m(m) { }

unsigned Pair::first() const { return n; }

unsigned Pair::second() const { return m; }

unsigned Pair::operator[](unsigned int index) const {
    if (index < 0 || index >= 2) throw "Invalid index";
    return (index == 0)? this->first(): this->second();
}
