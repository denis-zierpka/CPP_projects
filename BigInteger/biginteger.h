#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <complex>


using ld = long double;
using cld = std::complex<long double>;
long double pi = acos(-1);

std::vector<cld> convert_to_complex(const std::vector<int>& object, int max_size);
std::vector<int> convert_to_int(std::vector<cld>& object);
void fft(std::vector<cld>& object, bool reversed_fft);
void reversed_fft(std::vector<cld>& object);



class BigInteger {
public:
    static const int NUMBER_COUNT;
    static const int BASE;

private:
    std::vector<int> body_;
    bool negative_ = false;

public:
    BigInteger(): body_({0}), negative_(false) {}

    BigInteger(const BigInteger& other): body_(other.body_), negative_(other.negative_) {}

    BigInteger(const std::string& value) {
        body_.clear();
        if (value[0] == '-') {
            negative_ = true;
        }
        for (int i = static_cast<int>(value.size()) - 1; i >= (negative_ ? 1 : 0); --i) {
            int max_number_of_symbols = std::min(NUMBER_COUNT, i - (negative_ ? 1 : 0) + 1);
            int new_number = 0;
            int base_change = 1;
            for (int j = i; j > i - max_number_of_symbols; --j) {
                new_number += (value[j] - '0') * base_change;
                base_change *= 10;
            }
            i -= max_number_of_symbols - 1;
            body_.push_back(new_number);
        }
        normalize();
    }

    BigInteger(long long element): BigInteger(std::to_string(element)) {}

    BigInteger& operator= (const BigInteger& other) {
        if (this == &other) return *this;
        body_ = other.body_;
        negative_ = other.negative_;
        normalize();
        return *this;
    }

    std::string toString() const {
        std::string new_object;
        if (negative())
            new_object += '-';
        int max_size = size() - 1;
        for (int i = max_size; i >= 0; --i) {
            if (i != max_size) {
                std::string base_number = std::to_string(body_[i]);
                for (size_t j = 0; j < BigInteger::NUMBER_COUNT - base_number.size(); ++j)
                    new_object += '0';
            }
            new_object += std::to_string(body_[i]);
        }
        return new_object;
    }

    bool negative() const {
        return negative_;
    }

    void change_sign() {
        negative_ = !negative_;
    }

    int size() const {
        return static_cast<int>(body_.size());
    }

    BigInteger operator-() const {
        BigInteger new_object = *this;
        new_object.negative_ = !new_object.negative_;
        new_object.normalize();
        return new_object;
    }

    void normalize() {
        for (int i = 0; i < size(); ++i) {
            if (body_[i] >= BASE) {
                if (i + 1 < size()) {
                    body_[i + 1] += body_[i] / BASE;
                } else {
                    body_.push_back(body_[i] / BASE);
                }
                body_[i] %= BASE;
            }
            if (body_[i] < 0) {
                body_[i + 1] -= std::abs(body_[i]) / BASE;
                body_[i] %= BASE;
                if (body_[i] < 0) {
                    body_[i] += BASE;
                    body_[i + 1]--;
                }
            }
        }

        while (size() > 1 && body_[size() - 1] == 0)
            body_.pop_back();
        if (size() == 0) {
            body_.push_back(0);
            negative_ = false;
        } else if (size() == 1 && body_[0] == 0) {
            negative_ = false;
        }
    }

    BigInteger abs() const {
        BigInteger new_object = *this;
        if (new_object.negative()) {
            new_object.change_sign();
            new_object.normalize();
        }
        return new_object;
    }


    bool operator> (const BigInteger& other) const {
        if (negative_ == other.negative_) {
            return (negative_ ? other : *this).more_if_equal_sign(negative_ ? *this : other);
        }
        return !negative();
    }

    explicit operator bool() const {
        return *this != 0;
    }

    bool operator< (const BigInteger& other) const {
        return other > *this;
    }

    bool operator== (const BigInteger& other) const {
        return size() == other.size() && !(*this > other) && !(other > *this);
    }

    bool operator>= (const BigInteger& other) const {
        return *this > other || *this == other;
    }

    bool operator<= (const BigInteger& other) const {
        return *this < other || *this == other;
    }

    bool operator!= (const BigInteger& other) const {
        return !(*this == other);
    }

    BigInteger& operator+= (const BigInteger& other) {
        if (negative() == other.negative()) {
            plus_if_same_sign(other);
        } else {
            minus_if_same_sign(other);
        }
        return *this;
    }

    BigInteger& operator++() {
        *this += 1;
        normalize();
        return *this;
    }

    BigInteger operator++(int) {
        BigInteger new_object = *this;
        *this += 1;
        normalize();
        return new_object;
    }

    BigInteger& operator--() {
        *this -= 1;
        normalize();
        return *this;
    }

    BigInteger operator--(int) {
        BigInteger new_object = *this;
        *this -= 1;
        normalize();
        return new_object;
    }

    BigInteger& operator-= (const BigInteger& other) {
        if (negative() == other.negative()) {
            minus_if_same_sign(other);
        } else {
            plus_if_same_sign(other);
        }
        return *this;
    }

    BigInteger& operator*= (const BigInteger& other) {
        int max_size = std::max(size(), other.size());
        std::vector<cld> first_number = convert_to_complex(body_, max_size);
        std::vector<cld> second_number = convert_to_complex(other.body_, max_size);
        fft(first_number, false);
        fft(second_number, false);
        for (size_t i = 0; i < first_number.size(); ++i) {
            first_number[i] *= second_number[i];
        }
        reversed_fft(first_number);
        std::vector<int> new_number = convert_to_int(first_number);
        body_ = new_number;
        negative_ = (negative_ != other.negative_);
        normalize();
        return *this;
    }

    BigInteger& operator/= (const BigInteger& other) {
        BigInteger a = BigInteger();
        BigInteger new_body = BigInteger();
        division(other, a, new_body);
        *this = new_body;
        normalize();
        return *this;
    }

    BigInteger& operator%= (const BigInteger& other) {
        BigInteger a = BigInteger();
        BigInteger new_body = BigInteger();
        division(other, a, new_body);
        *this = a;
        normalize();
        return *this;
    }

private:
    bool more_if_equal_sign(const BigInteger& other) const {
        if (size() != other.size())
            return size() > other.size();

        for (int i = size() - 1; i >= 0; --i) {
            if (body_[i] < other.body_[i])
                return false;
            if (body_[i] > other.body_[i])
                return true;
        }
        return false;
    }

    void plus_if_same_sign(const BigInteger& other) {
        while (other.size() > size())
            body_.push_back(0);
        for (int i = 0; i < size() && i < other.size(); ++i) {
            body_[i] += other.body_[i];
        }
        normalize();
    }

    void minus_if_same_sign(const BigInteger& other) {
        bool second_more = false;
        if (!more_if_equal_sign(other)) {
            second_more = true;
        }
        while (other.size() > size())
            body_.push_back(0);
        for (int i = 0; i < size() && i < other.size(); ++i) {
            if (second_more) {
                body_[i] += other.body_[i] - 2 * body_[i];
            } else {
                body_[i] -= other.body_[i];
            }
        }
        if (second_more)
            negative_ = !negative_;
        normalize();
    }

    BigInteger& mult(const BigInteger& other) {
        int new_size = size() + other.size() + 2;
        BigInteger new_body = BigInteger();
        new_body.body_.resize(new_size);
        for (int i = 0; i < size(); ++i) {
            for (int j = 0; j < other.size(); ++j) {
                new_body.body_[i + j] += body_[i] * other.body_[j];
            }
        }
        new_body.normalize();
        new_body.negative_ = (negative_ != other.negative_);
        *this = new_body;
        normalize();
        return *this;
    }

    void division(const BigInteger& other, BigInteger& a, BigInteger& new_body) {
        if (*this == other) {
            a = 0;
            new_body = 1;
            return;
        }
        if (other == 0) {
            a = *this;
            new_body = *this;
            return;
        }
        if (other == 1 || other == -1) {
            a = 0;
            new_body = *this;
            new_body.negative_ = (negative_ != other.negative_);
            return;
        }
        a = 0;
        new_body = 0;
        for (int i = size() - 1; i >= 0; --i) {
            a.mult(BASE);
            a += body_[i];
            int q = 0;
            while (a.more_if_equal_sign(other) || a == other || a == -other) {
                if (!other.negative())
                    a -= other;
                else
                    a += other;
                ++q;
            }
            new_body.mult(BASE);
            new_body += q;
        }
        new_body.negative_ = (negative_ != other.negative_);
        a.negative_ = negative_;
    }
};

const int BigInteger::NUMBER_COUNT = 4;
const int BigInteger::BASE = 10000;

BigInteger operator"" _bi(unsigned long long element) {
    BigInteger new_object = element;
    return new_object;
}

BigInteger operator+ (const BigInteger& object, const BigInteger& other) {
    BigInteger new_object = object;
    new_object += other;
    return new_object;
}

BigInteger operator- (const BigInteger& object, const BigInteger& other) {
    BigInteger new_object = object;
    new_object -= other;
    return new_object;
}

BigInteger operator* (const BigInteger& object, const BigInteger& other) {
    BigInteger new_object = object;
    new_object *= other;
    return new_object;
}

BigInteger operator/ (const BigInteger& object, const BigInteger& other) {
    BigInteger new_object = object;
    new_object /= other;
    return new_object;
}

BigInteger operator% (const BigInteger& object, const BigInteger& other) {
    BigInteger new_object = object;
    new_object %= other;
    return new_object;
}

std::istream& operator >> (std::istream& in, BigInteger& input) {
    std::string new_input;
    in >> new_input;
    input = BigInteger(new_input);
    input.normalize();
    return in;
}

std::ostream& operator << (std::ostream& out, const BigInteger& output) {
    std::string result = output.toString();
    out << result;
    return out;
}

std::vector<cld> convert_to_complex(const std::vector<int>& object, int max_size) {
    int power_2 = 1;
    while (power_2 < max_size) {
        power_2 *= 2;
    }
    power_2 *= 2;
    std::vector<cld> ans;
    for (int i : object) {
        cld w(i);
        ans.push_back(w);
    }
    for (size_t i = 0; i < power_2 - object.size(); ++i) {
        cld w(0);
        ans.push_back(w);
    }
    return ans;
}

std::vector<int> convert_to_int(std::vector<cld>& object) {
    std::vector<int> ans;
    for (auto& i : object) {
        ans.push_back(i.real());
    }
    return ans;
}

int rb(int a, int n) {
    int ans = 0;
    while (a > 0) {
        ans += (a % 2) * (n / 2);
        a /= 2;
        n /= 2;
    }
    return ans;
}


void fft(std::vector<cld>& object, bool reversed_fft) {
    int n = object.size();
    for (int i = 0; i < n; ++i) {
        if (i < rb(i, n))
            swap(object[i], object[rb(i, n)]);
    }
    for (int len = 2; len <= n; len *= 2) {
        for (int i = 0; i < n; i += len) {
            ld angle = (reversed_fft ? -1 : 1) * 2 * pi / len;
            cld r = 1;
            cld w(cos(angle), sin(angle));
            for (int j = i; j < i + len / 2; ++j) {
                cld save = (object[j] + object[j + len / 2] * r);
                cld save2 = (object[j] - object[j + len / 2] * r);
                object[j] = save;
                object[j + len / 2] = save2;
                r *= w;
            }
        }
    }
}

void reversed_fft(std::vector<cld>& object) {
    fft(object, true);
    for (size_t i = 0; i < object.size(); ++i) {
        object[i] /= object.size();
        object[i] = floorl(object[i].real() + (ld)0.5);
    }
}



class Rational {
private:
    BigInteger numerator_;
    BigInteger denominator_;

public:
    Rational(): numerator_(0), denominator_(1) {}

    Rational(const Rational& other): numerator_(other.numerator_), denominator_(other.denominator_) {
        normalize();
    }

    Rational(const BigInteger& element): numerator_(element), denominator_(1) {}

    Rational(const int element): numerator_(element), denominator_(1) {}

    Rational& operator= (const Rational& other) {
        if (this == &other) return *this;
        numerator_ = other.numerator_;
        denominator_ = other.denominator_;
        normalize();
        return *this;
    }

    std::string toString() const {
        std::string new_object;
        new_object += numerator_.toString();
        if (denominator_ != 1) {
            new_object += '/';
            new_object += denominator_.toString();
        }
        return new_object;
    }

    std::string asDecimal(size_t precision = 0) const {
        BigInteger max_real = numerator_.abs() / denominator_;
        BigInteger rest = numerator_.abs() % denominator_;
        std::string new_object;
        if (numerator_.negative())
            new_object += '-';
        new_object += max_real.toString();
        if (precision <= 0)
            return new_object;

        new_object += '.';
        BigInteger base_number = 10;
        for (size_t i = 0; i < precision; ++i) {
            rest *= base_number;
            new_object += (rest / denominator_).toString();
            rest %= denominator_;
        }
        return new_object;
    }

    explicit operator double () const {
        std::stringstream stream;
        stream << asDecimal(100);
        double result;
        stream >> result;
        return result;
    }

    BigInteger gcd(const BigInteger& a, const BigInteger& b) {
        return a != 0 ? gcd(b % a, a) : b;
    }

    void normalize() {
        if (denominator_.negative()) {
            denominator_.change_sign();
            numerator_.change_sign();
        }
        numerator_.normalize();
        denominator_.normalize();
        BigInteger gcd_save = gcd(numerator_.abs(), denominator_.abs());
        if (gcd_save > 1) {
            numerator_ /= gcd_save;
            denominator_ /= gcd_save;
        }
    }

    Rational operator-() const {
        Rational new_object = *this;
        new_object.numerator_ = -new_object.numerator_;
        new_object.normalize();
        return new_object;
    }

    explicit operator bool() const {
        return *this != 0;
    }

    bool operator> (const Rational& other) const {
        return numerator_ * other.denominator_ > other.numerator_ * denominator_;
    }

    bool operator< (const Rational& other) const {
        return numerator_ * other.denominator_ < other.numerator_ * denominator_;
    }

    bool operator== (const Rational& other) const {
        return numerator_ == other.numerator_ && denominator_ == other.denominator_;
    }

    bool operator>= (const Rational& other) const {
        return *this > other || *this == other;
    }

    bool operator<= (const Rational& other) const {
        return *this < other || *this == other;
    }

    bool operator!= (const Rational& other) const {
        return !(*this == other);
    }

    Rational& operator+= (const Rational& other) {
        if (this == &other) {
            numerator_ *= 2;
            normalize();
            return *this;
        }
        numerator_ *= other.denominator_;
        numerator_ += other.numerator_ * denominator_;
        denominator_ *= other.denominator_;
        normalize();
        return *this;
    }

    Rational& operator-= (const Rational& other) {
        if (this == &other) {
            numerator_ = 0;
            normalize();
            return *this;
        }
        numerator_ *= other.denominator_;
        numerator_ -= other.numerator_ * denominator_;
        denominator_ *= other.denominator_;
        normalize();
        return *this;
    }

    Rational& operator*= (const Rational& other) {
        numerator_ *= other.numerator_;
        denominator_ *= other.denominator_;
        normalize();
        return *this;
    }

    Rational& operator/= (const Rational& other) {
        if (this == &other) {
            *this = 1;
            return *this;
        }
        numerator_ *= other.denominator_;
        denominator_ *= other.numerator_;
        normalize();
        return *this;
    }
};

Rational operator+ (const Rational& object, const Rational& other) {
    Rational new_object = object;
    new_object += other;
    return new_object;
}

Rational operator- (const Rational& object, const Rational& other) {
    Rational new_object = object;
    new_object -= other;
    return new_object;
}

Rational operator* (const Rational& object, const Rational& other) {
    Rational new_object = object;
    new_object *= other;
    return new_object;
}

Rational operator/ (const Rational& object, const Rational& other) {
    Rational new_object = object;
    new_object /= other;
    return new_object;
}