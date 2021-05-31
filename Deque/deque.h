#include <iostream>
#include <vector>

template <typename T>
class Deque {
private:
    std::vector<T> body;
    size_t L = 0;
    size_t R = 0;
    size_t sz = 0;
    int step_change = 0;

public:
    class iterator;
    class const_iterator;

    Deque<T>() = default;

    Deque<T>(size_t size, T object) {
        size_t new_size = (size % 2 == 0 ? size + 1 : size);
        L = 0;
        R = size;
        body.resize(new_size, object);
        sz = new_size;
    }

    Deque<T>(size_t size): Deque<T>(size, T()) {}

    Deque<T>& operator= (const Deque<T>& other) {
        body = other.body;
        L = other.L;
        R = other.R;
        sz = other.sz;
        return *this;
    }

    size_t size() const {
        return R - L;
    }

    T& operator[](size_t index) {
        return body[L + index];
    };

    T operator[](size_t index) const {
        return body[L + index];
    };

    T& at(size_t index) {
        if (index < 0 || index >= R - L)
            throw std::out_of_range("name");
        return body[L + index];
    };

    T at(size_t index) const {
        if (index < 0 || index >= R - L)
            throw std::out_of_range("name");
        return body[L + index];
    };

    void change_size() {
        std::vector<T> new_body(sz * 2 + 1);
        for (size_t i = 0; i < sz; ++i) {
            new_body[i + (sz - sz / 2)] = body[i];
        }
        L = L + (sz - sz / 2);
        R = R + (sz - sz / 2);
        sz = new_body.size();
        body = new_body;
    }

    void push_back(const T& object) {
        while (R >= sz)
            change_size();

        body[R] = object;
        ++R;
    }

    void push_front(const T& object) {
        while (L <= 0)
            change_size();
        --L;
        ++step_change;
        body[L] = object;
    }

    void pop_back() {
        --R;
    }

    void pop_front() {
        --step_change;
        ++L;
    }

    iterator begin() {
        iterator ans(this, 0, step_change);
        return ans;
    }

    iterator end() {
        iterator ans(this, R - L, step_change);
        return ans;
    }

    const_iterator begin() const {
        const_iterator ans(const_cast<Deque<T>*>(this), 0, step_change);
        return ans;
    }

    const_iterator end() const {
        const_iterator ans(const_cast<Deque<T>*>(this), R - L, step_change);
        return ans;
    }

    const_iterator cbegin() const {
        const_iterator ans(const_cast<Deque<T>*>(this), 0, step_change);
        return ans;
    }

    const_iterator cend() const {
        const_iterator ans(const_cast<Deque<T>*>(this), R - L, step_change);
        return ans;
    }

    void insert(const iterator& it, const T& object) {
        push_back(object);
        for (auto i = end() - 1; i > it; --i) {
            std::swap(*(i - 1), *i);
        }
    }
    void erase(const iterator& it) {
        for (auto i = it; i < end() - 1; ++i) {
            std::swap(*i, *(i + 1));
        }
        pop_back();
    }

    class iterator {
    public:
        Deque<T>* base;
        int it;
        int step;

        iterator() = default;

        iterator(Deque<T>* base, int it, int step): base(base), it(it), step(step) {}

        int real_position() const {
            return it - step;
        }

        iterator& operator++() {
            ++it;
            return *this;
        }

        iterator operator++(int) {
            iterator new_obj = *this;
            ++it;
            return new_obj;
        }

        iterator& operator--() {
            --it;
            return *this;
        }

        iterator operator--(int) {
            iterator new_obj = *this;
            --it;
            return new_obj;
        }

        iterator& operator+=(int value) {
            it += value;
            return *this;
        }

        iterator& operator-=(int value) {
            it -= value;
            return *this;
        }

        iterator operator+ (int value) const {
            iterator new_obj = *this;
            new_obj += value;
            return new_obj;
        }

        iterator operator- (int value) const {
            iterator new_obj = *this;
            new_obj -= value;
            return new_obj;
        }

        int operator-(const iterator& other) const {
            return real_position() - other.real_position();
        }

        bool operator< (const iterator& other) const {
            return real_position() < other.real_position();
        }

        bool operator> (const iterator& other) const {
            return real_position() > other.real_position();
        }

        bool operator== (const iterator& other) const {
            return real_position() == other.real_position();
        }

        bool operator<= (const iterator& other) const {
            return real_position() <= other.real_position();
        }

        bool operator>= (const iterator& other) const {
            return real_position() >= other.real_position();
        }

        bool operator!= (const iterator& other) const {
            return real_position() != other.real_position();
        }

        T& operator*() const {
            size_t L0 = base->L;
            int step_change0 = base->step_change;
            return base->body[L0 + step_change0 + real_position()];
        }

        T* operator->() const {
            size_t L0 = base->L;
            int step_change0 = base->step_change;
            return &base->body[L0 + step_change0 + real_position()];
        }
    };

    class const_iterator: public iterator {
    public:
        const_iterator(Deque<T>* base, int it, int step): iterator(base, it, step) {}

        const T& operator*() const {
            size_t L0 = iterator::base->L;
            int step_change0 = iterator::base->step_change;
            return iterator::base->body[L0 + step_change0 + iterator::real_position()];
        }

        const T* operator->() const {
            size_t L0 = iterator::base->L;
            int step_change0 = iterator::base->step_change;
            return &iterator::base->body[L0 + step_change0 + iterator::real_position()];
        }
    };
};
