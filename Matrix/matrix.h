#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <complex>
#include <cmath>
#include <bits/stl_algo.h>


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
    }

    BigInteger(long long element): BigInteger(std::to_string(element)) {}

    BigInteger& operator= (const BigInteger& other) {
        if (this == &other) return *this;
        body_ = other.body_;
        negative_ = other.negative_;
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
        if (new_object.size() == 1 && new_object.body_[0] == 0) {
            new_object.negative_ = false;
        }
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
        return *this;
    }

    BigInteger operator++(int) {
        BigInteger new_object = *this;
        *this += 1;
        return new_object;
    }

    BigInteger& operator--() {
        *this -= 1;
        return *this;
    }

    BigInteger operator--(int) {
        BigInteger new_object = *this;
        *this -= 1;
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
        if (size() < 100 || other.size() < 100) {
            mult(other);
            return *this;
        }
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
        if (abs() < other.abs()) {
            *this = BigInteger(0);
            negative_ = (negative_ != other.negative_);
            return *this;
        }
        division(other, a, new_body);
        *this = new_body;
        return *this;
    }

    BigInteger& operator%= (const BigInteger& other) {
        BigInteger a = BigInteger();
        BigInteger new_body = BigInteger();
        if (abs() < other.abs())
            return *this;
        division(other, a, new_body);
        *this = a;
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
            a.body_.push_back(0);
            for (int j = a.size() - 1; j > 0; --j)
                std::swap(a.body_[j], a.body_[j - 1]);

            a += body_[i];
            int q = 0;
            while (a.more_if_equal_sign(other) || a == other || a == -other) {
                if (!other.negative())
                    a -= other;
                else
                    a += other;
                ++q;
            }
            new_body.body_.push_back(0);
            for (int j = new_body.size() - 1; j > 0; --j)
                std::swap(new_body.body_[j], new_body.body_[j - 1]);
            new_body += q;
        }
        new_body.negative_ = (negative_ != other.negative_);
        a.negative_ = negative_;
    }

public:
    bool is_even() const {
        return body_[0] % 2 == 0;
    }

    void div_by_two() {
        std::vector<int> res;
        int qw = 0;
        for (int i = size() - 1; i >= 0; --i) {
            qw *= 10;
            qw += body_[i];
            if (!res.empty() || qw / 2 != 0)
                res.push_back(qw / 2);
            qw %= 2;
        }
        for (size_t i = 0; 2 * i < res.size(); ++i)
            std::swap(res[i], res[res.size() - i - 1]);

        body_ = res;
    }
};

const int BigInteger::NUMBER_COUNT = 1;
const int BigInteger::BASE = 10;

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

    Rational(const Rational& other): numerator_(other.numerator_), denominator_(other.denominator_) {}

    Rational(const BigInteger& element): numerator_(element), denominator_(1) {}

    Rational(const int element): numerator_(element), denominator_(1) {}

    Rational& operator= (const Rational& other) {
        if (this == &other) return *this;
        numerator_ = other.numerator_;
        denominator_ = other.denominator_;
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

    BigInteger gcd(BigInteger a, BigInteger b) {
        int two_divisors = 0;
        if (a == 0) {
            return b;
        } else if (b == 0) {
            return a;
        }
        BigInteger result = 1;
        while (a > 1 && b > 1) {
            if (a == b) {
                result = a;
                break;
            }
            if (a.is_even() && b.is_even()) {
                a.div_by_two();
                b.div_by_two();
                ++two_divisors;
            } else if (!a.is_even() && !b.is_even()) {
                if (a > b) {
                    a -= b;
                    a.div_by_two();
                } else {
                    b -= a;
                    b.div_by_two();
                }
            } else if (a.is_even()) {
                a.div_by_two();
            } else {
                b.div_by_two();
            }
        }
        for (int i = 0; i < two_divisors; ++i) {
            result *= 2;
        }
        return result;
    }

    void normalize() {
        if (denominator_.negative()) {
            denominator_.change_sign();
            numerator_.change_sign();
        }
        if (denominator_ == 1)
            return;
        BigInteger gcd_save = gcd(numerator_.abs(), denominator_.abs());
        if (gcd_save > 1) {
            numerator_ /= gcd_save;
            denominator_ /= gcd_save;
        }
    }

    Rational operator-() const {
        Rational new_object = *this;
        new_object.numerator_ = -new_object.numerator_;
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

    friend std::ostream& operator << (std::ostream& out, const Rational& output);
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

std::istream& operator >> (std::istream& in, Rational& input) {
    std::string new_input;
    in >> new_input;
    BigInteger d = BigInteger(new_input);
    input = Rational(d);
    input.normalize();
    return in;
}

std::ostream& operator << (std::ostream& out, const Rational& output) {
    out << output.numerator_ << '/' << output.denominator_;
    return out;
}



template<unsigned N, unsigned L, unsigned D>
struct square_helper {
    static const unsigned value = square_helper<N, (long long)(L + D / 2) * (L + D / 2) <= N ? (L + D / 2) : L,
            (long long)(L + D / 2) * (L + D / 2) <= N ? (D + 1) / 2 : D / 2>::value;
};

template<unsigned N, unsigned L>
struct square_helper<N, L, 1> {
    static const unsigned value = L;
};

template<unsigned N>
struct square {
    static const unsigned value = square_helper<N, 1, N>::value;
};

template<unsigned N>
static const unsigned square_v = square<N>::value;

template <unsigned N, unsigned M>
struct is_prime_helper {
    static const bool value = (N % M == 0 ? 0 : is_prime_helper<N, M - 1>::value);
};

template <unsigned N>
struct is_prime_helper<N, 1> {
    static const bool value = true;
};

template <unsigned N>
struct is_prime {
    static const bool value = is_prime_helper<N, square_v<N>>::value;
};

template <>
struct is_prime<0> {
    static const bool value = false;
};

template <>
struct is_prime<1> {
    static const bool value = false;
};

template <unsigned N>
static const bool is_prime_v = is_prime<N>::value;


template <unsigned N, unsigned M>
struct special_for {
    static const bool value = (((long long)N % (long long)M == 0 && (long long)M * (long long)M <= (long long)N) ?
                               (((long long)N / (long long)M) % (long long)M != 0 ? 0 : special_for<N, M - 1>::value) : special_for<N, M - 1>::value);
};

template <unsigned N>
struct special_for<N, 1> {
    static const bool value = true;
};

template <unsigned N>
struct is_prime_in_power {
    static const bool value = (N % 2 != 0 && special_for<N, square_v<N>>::value);
};

template <unsigned N>
struct has_primitive_root {
    static const bool value = is_prime_in_power<(N % 2 == 0) ? N / 2 : N>::value;
};

template <>
struct has_primitive_root<1> {
    static const bool value = false;
};

template <>
struct has_primitive_root<2> {
    static const bool value = true;
};

template <>
struct has_primitive_root<4> {
    static const bool value = true;
};

template <unsigned N>
static const bool has_primitive_root_v = has_primitive_root<N>::value;


template <unsigned N>
struct Exeption {};

template <>
struct Exeption<1> {
    static const bool value = true;
};

long long euler_function(unsigned N) {
    long long result = N, n = N;
    for (long long i = 2; i * i <= n; i++) {
        if (n % i == 0) {
            while (n % i == 0)
                n /= i;
            result -= result / i;
        }
    }
    if (n > 1) {
        result -= result / n;
    }
    return result;
}

template <unsigned N>
class Residue {
private:
    long long data;

public:
    Residue<N>() = default;

    Residue<N>(long long object): data(((object % (long long)N) + (long long)N) % (long long)N) {}

    explicit operator int() const {
        return data;
    }

    Residue<N> operator-() const {
        Residue<N> result = *this;
        result.data = (long long)(N - result.data) % (long long)N;
        return result;
    }

    Residue<N>& operator+= (const Residue<N>& other) {
        data += other.data;
        data %= (long long)N;
        return *this;
    }

    Residue<N>& operator-= (const Residue<N>& other) {
        data -= other.data;
        data = (data % (long long)N + (long long)N) % (long long)N;
        return *this;
    }

    Residue<N>& operator*= (const Residue<N>& other) {
        data *= other.data;
        data %= (long long)N;
        return *this;
    }

    bool operator==(const Residue<N>& other) const {
        return data == other.data;
    }

    explicit operator bool() const {
        return data != 0;
    }

    Residue<N>& operator=(const Residue<N>& other) {
        data = other.data;
        return *this;
    }

    bool operator!=(const Residue<N>& other) const {
        return !(*this == other);
    }

    Residue<N> pow(unsigned k) const {
        Residue<N> result(1);
        Residue<N> mult = *this;
        while (k) {
            if (k % 2 == 1)
                result *= mult;
            mult *= mult;
            k /= 2;
        }
        return result;
    }

    Residue<N> getInverse() const {
        static const unsigned b = is_prime_v<N>;
        Exeption<b> a;
        if (a.value) {
            Residue<N> result = *this;
            result = result.pow(N - 2);
            return result;
        }
    }

    Residue<N> operator/ (const Residue<N>& other) const {
        static const unsigned b = is_prime_v<N>;
        Exeption<b> a;
        if (a.value) {
            Residue<N> result = *this;
            result *= other.getInverse();
            return result;
        }
    }


    long long order() const {
        long long phi = euler_function(N);
        long long last = 0;
        for (long long i = 1; i * i <= phi; i++) {
            if (phi % i == 0 && pow(i).data == 1) {
                return i;
            }
            if (phi % i == 0 && pow(phi / i).data == 1) {
                last = phi / i;
            }
        }
        return last;
    }

    static Residue<N> getPrimitiveRoot() {
        long long phi = euler_function(N);
        for (unsigned res = 2; res <= N; ++res) {
            if (std::__gcd((long long)res, (long long)N) != 1 || Residue<N>(res).pow(phi / 2).data == 1)
                continue;
            bool ok = true;
            Residue<N> re(res);
            for (unsigned i = 2; i * i <= phi; ++i) {
                if (phi % i == 0) {
                    if (re.pow(i).data == 1 || re.pow(phi / i).data == 1) {
                        ok = false;
                        break;
                    }
                }
            }
            if (ok) {
                return re;
            }
        }
        Residue<N> result(0);
        return result;
    }
};


template <unsigned N>
Residue<N> operator+ (const Residue<N>& first, const Residue<N>& second) {
    Residue<N> result = first;
    result += second;
    return result;
}

template <unsigned N>
Residue<N> operator- (const Residue<N>& first, const Residue<N>& second) {
    Residue<N> result = first;
    result -= second;
    return result;
}

template <unsigned N>
Residue<N> operator* (const Residue<N>& first, const Residue<N>& second) {
    Residue<N> result = first;
    result *= second;
    return result;
}


template <unsigned M, unsigned N, typename Field = Rational>
class Matrix {
public:
    std::vector<std::vector<Field>> body;

    void arithmetic_operation(const Matrix<M, N, Field> other, bool subtraction) {
        for (unsigned i = 0; i < M; ++i) {
            for (unsigned j = 0; j < N; ++j) {
                body[i][j] += Field(subtraction ? -1 : 1) * other.body[i][j];
            }
        }
    }

    void sum_of_rows(int a, int b, int start, Field coff) {
        for (unsigned j = start; j < N; ++j) {
            body[b][j] += body[a][j] * coff;
        }
    }

    void stepped_matrix_view() {
        int index = 0;
        for (unsigned i = 0; i < N; ++i) {
            int find = -1;
            for (unsigned j = index; j < M; ++j) {
                if (body[j][i] != 0) {
                    find = j;
                    break;
                }
            }
            if (find == -1)
                continue;
            std::swap(body[index], body[find]);
            for (unsigned j = index + 1; j < M; ++j) {
                sum_of_rows(index, j, 0, -body[j][i] / body[index][i]);
            }
            index++;
        }
    }

    template<unsigned K>
    Matrix<M, N + K, Field> concatenation_right(const Matrix<M, K, Field>& other) {
        Matrix<M, N + K, Field> save;
        save.body.clear();
        save.body.resize(M);
        for (unsigned i = 0; i < M; ++i) {
            for (unsigned j = 0; j < N; ++j) {
                save.body[i].push_back(body[i][j]);
            }
        }
        for (unsigned i = 0; i < M; ++i) {
            for (unsigned j = 0; j < K; ++j) {
                save.body[i].push_back(other.body[i][j]);
            }
        }
        return save;
    }

    void divide_the_row(int index, Field object) {
        for (unsigned j = index; j < N; ++j) {
            body[index][j] = body[index][j] / object;
        }
    }

    Matrix<M, N - M, Field> delete_left() const {
        Matrix<M, N - M, Field> res;
        res.body.clear();
        res.body.resize(M);
        for (unsigned i = 0; i < M; ++i) {
            for (unsigned j = M; j < N; ++j) {
                res.body[i].push_back(body[i][j]);
            }
        }
        return res;
    }

    Matrix<M, N, Field>() {
        body.resize(M, std::vector<Field>(N));
    }

    Matrix<M, N, Field> (std::initializer_list<std::initializer_list<Field>> list) {
        body.clear();
        body.resize(M);
        int i = 0;
        for (std::initializer_list<Field> i2 : list) {
            for (Field j2 : i2) {
                body[i].push_back(j2);
            }
            ++i;
        }
    }

    Matrix<M, N, Field>& operator= (const Matrix<M, N, Field>& other) {
        for (unsigned i = 0; i < M; ++i) {
            for (unsigned j = 0; j < N; ++j) {
                body[i][j] = other.body[i][j];
            }
        }
        return *this;
    }

    bool operator== (const Matrix<M, N, Field>& other) const {
        for (unsigned i = 0; i < M; ++i) {
            for (unsigned j = 0; j < N; ++j) {
                if (body[i][j] != other.body[i][j])
                    return false;
            }
        }
        return true;
    }

    bool operator!= (const Matrix<M, N, Field>& other) const {
        return !(*this == other);
    }

    Matrix<M, N, Field>& operator+= (const Matrix<M, N, Field>& other) {
        arithmetic_operation(other, false);
        return *this;
    }

    Matrix<M, N, Field>& operator-= (const Matrix<M, N, Field>& other) {
        arithmetic_operation(other, true);
        return *this;
    }

    Matrix<M, N, Field>& operator*= (Field object) {
        for (unsigned i = 0; i < M; ++i) {
            for (unsigned j = 0; j < N; ++j) {
                body[i][j] *= object;
            }
        }
        return *this;
    }

    Matrix<M, N, Field> operator+ (const Matrix<M, N, Field>& other) const {
        Matrix<M, N, Field> res = *this;
        res += other;
        return res;
    }

    Matrix<M, N, Field> operator- (const Matrix<M, N, Field>& other) const {
        Matrix<M, N, Field> res = *this;
        res -= other;
        return res;
    }

    Matrix<M, N, Field> operator* (Field object) const {
        Matrix<M, N, Field> res = *this;
        res *= object;
        return res;
    }

    Matrix<M, M, Field>& operator*= (const Matrix<M, M, Field>& other) {
        Matrix<M, M, Field> result = *this * other;
        *this = result;
        return *this;
    }

    template<unsigned K>
    Matrix<M, K, Field> operator* (const Matrix<N, K, Field>& other) const {
        Matrix<M, K, Field> result;
        result.body.resize(M, std::vector<Field>(K));
        for (unsigned m = 0; m < M; ++m) {
            for (unsigned k = 0; k < K; ++k) {
                for (unsigned n = 0; n < N; ++n) {
                    result.body[m][k] += body[m][n] * other.body[n][k];
                }
            }
        }
        return result;
    }

    Matrix<N, M, Field> transposed() const {
        Matrix<N, M, Field> res;
        res.body.resize(N, std::vector<Field>(M));
        for (unsigned m = 0; m < M; ++m) {
            for (unsigned n = 0; n < N; ++n) {
                res.body[n][m] = body[m][n];
            }
        }
        return res;
    }


    std::vector<Field> getRow(unsigned row) const {
        std::vector<Field> result;
        for (unsigned n = 0; n < N; ++n) {
            result.push_back(body[row][n]);
        }
        return result;
    }

    std::vector<Field> getColumn(unsigned column) const {
        std::vector<Field> result;
        for (unsigned m = 0; m < M; ++m) {
            result.push_back(body[m][column]);
        }
        return result;
    }

    std::vector<Field>& operator[] (size_t index) {
        return body[index];
    }

    std::vector<Field> operator[] (size_t index) const {
        return body[index];
    }

    Field det() const {
        static_assert(N == M);
        Matrix<M, N, Field> check = *this;
        check.stepped_matrix_view();
        Field result = 1;
        for (unsigned i = 0; i < N; ++i) {
            result *= check.body[i][i];
        }
        return result;
    }

    int rank() const {
        Matrix<M, N, Field> check = *this;
        check.stepped_matrix_view();
        int cnt = 0;
        for (unsigned i = 0; i < M; ++i) {
            int fl = 0;
            for (unsigned j = 0; j < N; ++j) {
                if (check.body[i][j] != 0) {
                    fl = 1;
                    break;
                }
            }
            if (fl == 1)
                cnt++;
        }

        return cnt;
    }

    Matrix<M, N, Field> inverted() const {
        static_assert(N == M);
        Matrix<M, N, Field> result = *this;
        result.invert();
        return result;
    }

    Matrix<M, N, Field>& invert() {
        static_assert(N == M);
        Matrix<M, M, Field> qwe;
        qwe.body.resize(M, std::vector<Field>(M));
        for (unsigned i = 0; i < M; ++i)
            qwe.body[i][i] = 1;

        Matrix<M, N + M, Field> first = concatenation_right(qwe);
        first.stepped_matrix_view();

        for (unsigned i = 0; i < M; ++i) {
            if (first.body[i][i] != 0)
                first.divide_the_row(i, first.body[i][i]);
        }

        for (unsigned j = (int)M - 1; j > 0; --j) {
            if (first.body[j][j] != 0) {
                for (unsigned k = 0; k < j; ++k) {
                    if (first.body[k][j] != 0) {
                        first.sum_of_rows(j, k, N, -first.body[k][j] / first.body[j][j]);
                        first.body[k][j] = 0;
                    }
                }
            }
        }
        *this = first.delete_left();
        return *this;
    }

    Field trace() const {
        static_assert(N == M);
        Field res = 0;
        for (unsigned i = 0; i < N; ++i) {
            res += body[i][i];
        }
        return res;
    }
};

template <unsigned M, unsigned N, typename Field = Rational>
Matrix<M, N, Field> operator* (Field object, const Matrix<M, N, Field>& other) {
    Matrix<M, N, Field> res = other;
    res *= object;
    return res;
}

template<unsigned N, typename Field = Rational>
using SquareMatrix = Matrix<N, N, Field>;
