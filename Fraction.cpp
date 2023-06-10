//
// Created by squir on 19.03.2023.
//

#include "Fraction.h"
Fraction abs(Fraction &f) {
    return f < 0 ? -f : f;
}

std::ostream &operator<<(std::ostream &out, const Fraction &f) {
    if (f.enumerator() == 0)
        out << 0;
    else {
        out << f.enumerator();
        if (f.denominator() != 1)
            out << '/' << f.denominator();
    }
    return out;
}

std::istream &operator>>(std::istream &in, Fraction &f) {
    std::string s;
    in >> s;
    auto pos = s.find('/');
    if (pos != std::string::npos) {
        std::string enumerator(s.begin(), std::next(s.begin(),
                                                    pos));
        std::string denominator(std::next(s.begin(), pos + 1),
                                s.end());
        f = {std::stoi(enumerator), std::stoi(denominator)};
    } else {
        f = std::stoi(s);
    }
    return in;
}
