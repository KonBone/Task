////
//// Created by konst on 06.09.2022.
////
//
//#include "Complex.h"
//
//Complex::Complex(double real, double imaginary) {
//    this->real = real;
//    this->imaginary = imaginary;
//}
//
//Complex & Complex::operator+(Complex* that) {
//    double re = this->real + that->real;
//    double im = this->imaginary + that->imaginary;
//    Complex *res = new Complex(re, im);
//    return *res;
//}

#include "Complex.h"

Complex::Complex(double real) {

}

Complex::Complex(double real, double imaginary) {

}

Complex::Complex(Complex &that) {

}

double Complex::getReal() {
    return 0;
}

double Complex::getImaginary() {
    return 0;
}

double Complex::getArgument() {
    return 0;
}

double Complex::abs() {
    return 0;
}

std::string Complex::toString() {
    return std::string();
}

Complex Complex::operator+(Complex that) {
    return 0;
}

Complex Complex::operator+(void) {
    return 0;
}

void Complex::operator+=(Complex that) {

}

Complex Complex::operator-(Complex that) {
    return 0;
}

Complex Complex::operator-(void) {
    return 0;
}

void Complex::operator-=(Complex that) {

}

Complex Complex::operator*(Complex that) {
    return 0;
}

void Complex::operator*=(Complex that) {

}

Complex Complex::operator/(Complex that) {
    return 0;
}

void Complex::operator/=(Complex that) {

}

//Complex &Complex::operator=(Complex that) {
//    return 0;
//}

