//
// Created by konst on 06.09.2022.
//

#ifndef TASK_COMPLEX_H
#define TASK_COMPLEX_H

#include <string>

class Complex {
private:
    double real, imaginary;
public:
    Complex(double real);
    Complex(double real, double imaginary);
    Complex(Complex & that);

    double getReal();
    double getImaginary();

    double getArgument();
    double abs();

    std::string toString();

    Complex operator +(Complex that);
    Complex operator +(void);
    void operator += (Complex that);

    Complex operator -(Complex that);
    Complex operator -(void);
    void operator -= (Complex that);

    Complex operator *(Complex that);
    void operator *= (Complex that);

    Complex operator /(Complex that);
    void operator /= (Complex that);

    Complex & operator = (Complex that);
};


#endif //TASK_COMPLEX_H
