#include <iostream>
#include <cstring>

class String;
bool operator==(const String& a, const String& b);

class String {
private:
    size_t capacity_;
    size_t length_;
    char* body_;

public:
    String(): capacity_(1), length_(0), body_(new char [1]) {}

    String(size_t n, char element): capacity_(n),
            length_(n), body_(new char [capacity_]) {

        memset(body_, element, length_);
    }

    String(const String& other): capacity_(other.capacity_),
            length_(other.length_), body_(new char [other.capacity_]) {

        memcpy(body_, other.body_, length_);
    }

    String(const char *list): capacity_(strlen(list)),
                              length_(strlen(list)), body_(new char [capacity_]) {

        memcpy(body_, list, length_);
    }

    String(char element): capacity_(1), length_(1),
                                body_(new char [1]) {
        body_[0] = element;
    }

    ~String() {
        delete[] body_;
    }

    void push_back(char element) {
        if (length_ == capacity_) {
            char* new_body = new char [capacity_ * 2];
            memcpy(new_body, body_, length_);
            capacity_ *= 2;
            delete[] body_;
            body_ = new_body;
        }
        body_[length()] = element;
        length_++;
    }

    void pop_back() {
        --length_;
    }

    size_t length() const {
        return length_;
    }

    char& front() {
        return body_[0];
    }

    char front() const {
        return body_[0];
    }

    char& back() {
        return body_[length_ - 1];
    }

    char back() const {
        return body_[length_ - 1];
    }

    bool empty() const {
        return (length() == 0);
    }

    void clear() {
        length_ = 0;
    }

    void swap(String& other) {
        std::swap(body_, other.body_);
        std::swap(length_, other.length_);
        std::swap(capacity_, other.capacity_);
    }

    String& operator= (const String& other) {
        if (this == &other) return *this;
        String new_other = other;
        swap(new_other);
        return *this;
    }

    char& operator[] (size_t index) {
        return body_[index];
    }

    char operator[] (size_t index) const {
        return body_[index];
    }

    String& operator+= (char element) {
        push_back(element);
        return *this;
    }

    String& operator+= (const String& other) {
        for (size_t i = 0; i < other.length(); ++i)
            push_back(other[i]);
        return *this;
    }

    bool equal(const String& other) const {
        if (length() != other.length()) return false;
        for (size_t i = 0; i < length(); ++i) {
            if (body_[i] != other[i])
                return false;
        }
        return true;
    }

    String substr(size_t start, size_t count) const {
        String return_value = String();
        return_value.capacity_ = count;
        delete[] return_value.body_;
        return_value.body_ = new char [return_value.capacity_];
        return_value.length_ = count;
        memcpy(return_value.body_, (body_ + start), count);
        return return_value;
    }

    size_t find(const String& substring) const {
        for (size_t i = 0; i < length() - substring.length() + 1; ++i) {
            if (substr(i, substring.length()) == substring) {
                return i;
            }
        }
        return length();
    }

    size_t rfind(const String& substring) const {
        for (size_t i = length() - substring.length();; --i) {
            if (substr(i, substring.length()) == substring) {
                return i;
            }
            if (i == 0)
                break;
        }
        return length();
    }
};

bool operator==(const String& a, const String& b) {
    return a.equal(b);
}

String operator+ (const String& a, const String& b) {
    String return_value = a;
    return_value += b;
    return return_value;
}

std::istream& operator >> (std::istream &cin, String& a) {
    char ch;
    while (cin >> std::noskipws >> ch) {
        if (ch != ' ' && ch != '\n')
            a.push_back(ch);
        else
            return cin;
    }
    return cin;
}

std::ostream& operator << (std::ostream &cout, const String& a) {
    for (size_t i = 0; i < a.length(); ++i)
        cout << a[i];
    return cout;
}
