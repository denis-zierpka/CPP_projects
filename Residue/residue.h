#include <iostream>
#include <cmath>
#include <bits/stl_algo.h>


template<unsigned N, unsigned L, unsigned D>
struct square_helper {
    static const unsigned value = square_helper<N, static_cast<long long>(L + D / 2) * (L + D / 2) <= N ? (L + D / 2) : L,
            static_cast<long long>(L + D / 2) * (L + D / 2) <= N ? (D + 1) / 2 : D / 2>::value;
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
    static const bool value = ((N % static_cast<long long>(M) == 0 && static_cast<long long>(M) * M <= N) ?
            ((static_cast<long long>(N) / M) % static_cast<long long>(M) != 0 ? 0 : special_for<N, M - 1>::value) : special_for<N, M - 1>::value);
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
struct Exception {};

template <>
struct Exception<1> {
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

    explicit Residue<N>(long long object): data(((object % static_cast<long long>(N)) + N) % static_cast<long long>(N)) {}

    explicit operator int() const {
        return data;
    }

    Residue<N>& operator+= (const Residue<N>& other) {
        data += other.data;
        data %= static_cast<long long>(N);
        return *this;
    }

    Residue<N>& operator-= (const Residue<N>& other) {
        data -= other.data;
        data = (data % static_cast<long long>(N) + N) % static_cast<long long>(N);
        return *this;
    }

    Residue<N>& operator*= (const Residue<N>& other) {
        data *= other.data;
        data %= static_cast<long long>(N);
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
        Exception<b> a;
        if (a.value) {
            Residue<N> result = *this;
            result = result.pow(N - 2);
            return result;
        }
    }

    Residue<N> operator/ (const Residue<N>& other) const {
        static const unsigned b = is_prime_v<N>;
        Exception<b> a;
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
            if (std::__gcd(static_cast<long long>(res), static_cast<long long>(N)) != 1 || Residue<N>(res).pow(phi / 2).data == 1)
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