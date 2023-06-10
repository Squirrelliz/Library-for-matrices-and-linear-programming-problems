//
// Created by squir on 19.03.2023.
//

#ifndef IOP_FRACTION_H
#define IOP_FRACTION_H

#include <iostream>
#include <ostream>
#include <numeric>

//класс Дробь
class Fraction {
public:
    Fraction() = default;

    //конструктор с входными данными числитель enumerator и знаменатель denominator
    Fraction(int enumerator, int denominator) {
        if (denominator == 0)
            throw std::domain_error("Division by zero cannot be performed\n");

        _enumerator = enumerator;
        _denominator = denominator;
        calculate();
    }

    //конструктор целого числа в формате дроби со значением enumerator
    Fraction(const int enumerator) {
        _enumerator = enumerator;
        _denominator = 1;
    }

    explicit operator int() const { return _enumerator / _denominator; }

    Fraction operator+(const Fraction &f) const {
        Fraction res;
        res._enumerator = _enumerator * f._denominator +
                          f._enumerator * _denominator;
        res._denominator = _denominator * f._denominator;

        res.calculate();

        return res;
    }

    Fraction &operator+=(const Fraction &f) {
        *this = *this + f;
        return *this;
    }

    Fraction operator+(const int number) const {
        Fraction res;
        res._enumerator = _enumerator + number * _denominator;
        res._denominator = _denominator;

        res.calculate();

        return res;
    }



    Fraction &operator+=(const int number) {
        _enumerator += number * _denominator;
        calculate();

        return *this;
    }

    Fraction operator-(const Fraction &f) const {
        Fraction res;
        res._enumerator = _enumerator * f._denominator -
                          f._enumerator * _denominator;
        res._denominator = _denominator * f._denominator;
        res.calculate();

        return res;
    }

    Fraction &operator-=(const Fraction &f) {
        *this = *this - f;
        return *this;
    }

    Fraction operator-(int number) const {
        Fraction res;
        res._enumerator = _enumerator - number * _denominator;
        res._denominator = _denominator;

        res.calculate();

        return res;
    }

    Fraction &operator-=(int num) {
        _enumerator -= num * _denominator;
        calculate();
        return *this;
    }

    Fraction operator*(const Fraction &f) const {
        Fraction res;
        res._enumerator = _enumerator * f._enumerator;
        res._denominator = _denominator * f._denominator;

        res.calculate();

        return res;
    }

    Fraction &operator*=(Fraction &f) {
        *this = *this * f;
        return *this;
    }

    Fraction operator*(int number) const {
        Fraction res;
        res._enumerator = _enumerator * number;
        res._denominator = _denominator;

        res.calculate();

        return res;
    }

    Fraction &operator*=(int number) {
        _enumerator *= number;
        calculate();
        return *this;
    }

    Fraction operator/(Fraction &f) const {
        Fraction res;
        res._enumerator = _enumerator * f._denominator;
        res._denominator = _denominator * f._enumerator;

        res.calculate();

        return res;
    }

    Fraction &operator/=(Fraction &f) {
        *this = *this / f;
        return *this;
    }

    Fraction operator/(int number) const {
        Fraction res;
        res._enumerator = _enumerator;
        res._denominator = _denominator * number;

        res.calculate();

        return res;
    }

    Fraction &operator/=(int number) {
        _denominator *= number;
        calculate();
        return *this;
    }

    Fraction &operator=(int number) {
        _enumerator = number;
        _denominator = 1;
        return *this;
    }

    Fraction operator-() {
        Fraction f = *this;
        f._enumerator = -_enumerator;
        return f;
    }

    bool operator==(const Fraction &f) const {
        return _enumerator == f._enumerator && _denominator == f._denominator;
    }

    bool operator==(const int number) const {
        return _enumerator == number;
    }

    bool operator!=(const Fraction &f) const {
        return _enumerator != f._enumerator || _denominator != f._denominator;
    }

    bool operator!=(const int number) const {
        return _enumerator != number;
    }

    bool operator>(const Fraction &f) const {
        int enumerator1 = _enumerator * f._denominator;
        int enumerator2 = f._enumerator * _denominator;
        return enumerator1 > enumerator2;
    }

    bool operator<(const Fraction &f) const {
        int enumerator1 = _enumerator * f._denominator;
        int enumerator2 = f._enumerator * _denominator;
        return enumerator1 < enumerator2;
    }

    bool operator>=(const Fraction &f) const {
        int enumerator1 = _enumerator * f._denominator;
        int enumerator2 = f._enumerator * _denominator;
        return enumerator1 >= enumerator2;
    }

    bool operator<=(const Fraction &f) const {
        int enumerator1 = _enumerator * f._denominator;
        int enumerator2 = f._enumerator * _denominator;
        return enumerator1 <= enumerator2;
    }


    bool operator>(const int number) const {
        return _enumerator > number * _denominator;
    }

    bool operator<(const int num) const {
        return _enumerator < num * _denominator;
    }

    bool operator>=(const int number) const {
        return _enumerator >= number * _denominator;
    }

    bool operator<=(const int number) const {
        return _enumerator <= number * _denominator;
    }

    //возвращает числитель дроби enumerator()
    [[nodiscard]] int enumerator() const {
        return _enumerator;
    }

    //возвращает знаменатель дроби denominator()
    [[nodiscard]] int denominator() const {
        return _denominator;
    }

private:
    int _enumerator = 0;//числитель
    int _denominator = 1;//знаменатель

    //выполняет сокращение дроби
    void calculate() {
        int gcd = std::gcd(_enumerator, _denominator);
        _enumerator /= gcd;
        _denominator /= gcd;
        if (_enumerator < 0 && _denominator < 0 ||
            _enumerator >= 0 && _denominator < 0) {
            _enumerator = -_enumerator;
            _denominator = -_denominator;
        }
    }
};

//возвращает модуль дроби f()
Fraction abs(Fraction &f);

std::ostream &operator<<(std::ostream &out, const Fraction &f);

std::istream &operator>>(std::istream &in, Fraction &f);

#endif //IOP_FRACTION_H
