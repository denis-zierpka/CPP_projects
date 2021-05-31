#include <iostream>
#include <vector>
#include <cmath>

long double PI = acos(-1);
long double EPS = 0.00001;

struct Point {
    long double x, y;

    Point(): x(0), y(0) {}

    Point(long double x, long double y): x(x), y(y) {}

    bool operator==(const Point& other) const {
        return (std::abs(x - other.x) < EPS && std::abs(y - other.y) < EPS);
    }
    bool operator!=(const Point& other) const {
        return !(*this == other);
    }
};


struct Vector {
    long double x, y;

    Vector(): x(0), y(0) {}

    Vector(const Point& a, const Point& b): x(b.x - a.x), y(b.y - a.y) {}

    Vector(long double a, long double b): x(a), y(b) {}

    long double size() const {
        return sqrtl(x * x + y * y);
    }

    Vector operator*= (long double value) {
        x *= value;
        y *= value;
        return *this;
    }

    Vector operator/= (long double value) {
        x /= value;
        y /= value;
        return *this;
    }

    void rotate(long double angle) {
        long double si = sinl(PI * angle / 180);
        long double co = cosl(PI * angle / 180);
        long double new_x = co * x - si * y;
        long double new_y = si * x + co * y;
        x = new_x;
        y = new_y;
    }
};

Point operator+ (const Point& a, const Vector& b) {
    Point new_point = Point(a.x + b.x, a.y + b.y);
    return new_point;
}

Point middle(const Point& a, const Point& b) {
    Vector v(a, b);
    v /= 2;
    Point ans = a + v;
    return ans;
}

long double scalar_product(const Vector& a, const Vector& b) {
    return a.x * b.x + a.y * b.y;
}

long double vector_product(const Vector& a, const Vector& b) {
    return a.x * b.y - a.y * b.x;
}

long double space(const Point& a, const Point& b) {
    Vector c = Vector(a, b);
    return c.size();
}

class Line {
private:
    long double A;
    long double B;
    long double C;

public:
    Line(const Point& first, const Point& second): A(second.y - first.y), B(first.x - second.x),
            C(-A * first.x - B * first.y) {}

    Line(long double k, long double b): A(k), B(-1), C(b) {}

    Line(const Point& first, long double k): Line(k, first.y - k * first.x) {}

    bool operator== (const Line& other) const {
        return std::abs(A * other.B - B * other.A) < EPS && std::abs(A * other.C - C * other.A) < EPS && std::abs(B * other.C - C * other.B) < EPS;
    }

    bool operator!= (const Line& other) const {
        return !(*this == other);
    }

    long double hight(const Point& object) const {
        return std::abs(A * object.x + B * object.y + C) / sqrtl(A * A + B * B);
    }

    Vector normal_vector() const {
        Vector v(A, B);
        return v;
    }

    long double result(const Point& object) const {
        return A * object.x + B * object.y + C;
    }

    Point intersection_point(const Line& l2) const {
        return Point((B * l2.C - l2.B * C) / (A * l2.B - l2.A * B), (l2.A * C - A * l2.C) / (A * l2.B - l2.A * B));
    }
};

class Shape {
public:
    virtual ~Shape() = default;
    virtual long double perimeter() const = 0;
    virtual long double area() const = 0;
    virtual bool operator== (const Shape& another) const = 0;
    virtual bool operator!= (const Shape& another) const = 0;
    virtual bool isCongruentTo(const Shape& another) const = 0;
    virtual bool isSimilarTo(const Shape& another) const = 0;
    virtual bool containsPoint(const Point& point) const = 0;

    virtual void rotate(const Point& center, long double angle) = 0;
    virtual void reflex(const Point& center) = 0;
    virtual void reflex(const Line& axis) = 0;
    virtual void scale(const Point& center, long double coefficient) = 0;
};

class Polygon: public Shape {
protected:
    std::vector<Point> Vertices;

public:
    Polygon(std::vector<Point>& object): Vertices(object) {}

    Polygon(std::initializer_list<Point> list): Vertices(list) {}

    int verticesCount() const {
        return Vertices.size();
    }
    std::vector<Point> getVertices() const {
        return Vertices;
    }

    bool isConvex() const {
        bool all_neg = true;
        bool all_pos = true;
        for (int i = 0; i < verticesCount(); ++i) {
            Vector edge1 = Vector(Vertices[i], Vertices[(i + 1) % verticesCount()]);
            Vector edge2 = Vector(Vertices[(i + 1) % verticesCount()], Vertices[(i + 2) % verticesCount()]);
            long double result = vector_product(edge1, edge2);
            if (result > 0) {
                all_neg = false;
            } else if (result < 0) {
                all_pos = false;
            }
        }
        return (all_pos != all_neg);
    }

    void add_Point(const Point& a) {
        Vertices.push_back(a);
    }

    long double perimeter() const override  {
        long double perimeter_length = 0;
        for (int i = 0; i < verticesCount(); ++i) {
            Vector edge = Vector(Vertices[i], Vertices[(i + 1) % verticesCount()]);
            perimeter_length += edge.size();
        }
        return perimeter_length;
    }

    long double area() const override {
        long double area_size = 0;
        for (int i = 1; i < verticesCount(); ++i) {
            Vector edge1 = Vector(Vertices[0], Vertices[i]);
            Vector edge2 = Vector(Vertices[0], Vertices[(i + 1) % verticesCount()]);
            area_size += vector_product(edge1, edge2);
        }
        area_size /= (long double)2;
        return area_size;
    }

    bool operator== (const Shape& another) const override {
        try {
            const auto& other = dynamic_cast<const Polygon&>(another);

            if (verticesCount() != other.verticesCount()) {
                return false;
            }
            Polygon new_this = *this;
            for (int ind = 0; ind < other.verticesCount(); ++ind) {
                bool fl = true;
                for (int i = 0; i < new_this.verticesCount(); ++i) {
                    if (new_this.Vertices[i] != other.Vertices[(i + ind) % other.verticesCount()]) {
                        fl = false;
                        break;
                    }
                }
                if (fl)
                    return true;
            }
            for (int ind = 0; ind < other.verticesCount(); ++ind) {
                bool fl = true;
                for (int i = 0; i < new_this.verticesCount(); ++i) {
                    if (new_this.Vertices[new_this.verticesCount() - i - 1] != other.Vertices[(i + ind) % other.verticesCount()]) {
                        fl = false;
                        break;
                    }
                }
                if (fl)
                    return true;
            }
            return false;
        } catch (std::bad_cast& b) {
            return false;
        }
        return false;
    }

    bool operator!= (const Shape& another) const override {
        return !(*this == another);
    }

private:
    std::vector<Point> delete_on_one_line(std::vector<Point> old) const {
        std::vector<Point> new_this;
        new_this.push_back(old[0]);
        for (size_t i = 1; i < old.size(); ++i) {
            if (new_this.size() < 2) {
                new_this.push_back(old[i]);
                continue;
            }
            Vector v = Vector(new_this[new_this.size() - 2], new_this[new_this.size() - 1]);
            Vector v2 = Vector(new_this[new_this.size() - 1], old[i]);
            if (std::abs(vector_product(v, v2) / v.size() / v2.size()) < EPS) {
                new_this.pop_back();
            }
            new_this.push_back(old[i]);
        }
        return new_this;
    }

    bool check_congruency(const Polygon& new_this, const Polygon& new_other, bool reversed_order) const {
        int save = new_other.verticesCount() - 1;
        for (int ind = 0; ind < new_other.verticesCount(); ++ind) {
            bool fl = true;
            for (int i = 0; i < new_this.verticesCount(); ++i) {
                Point this_1 = new_this.Vertices[i];
                Point this_2 = new_this.Vertices[(i + 1) % new_this.verticesCount()];
                Point this_3 = new_this.Vertices[(i + 2) % new_this.verticesCount()];

                Point other_1 = new_other.Vertices[(i + ind) % new_other.verticesCount()];
                Point other_2 = new_other.Vertices[(i + ind + 1) % new_other.verticesCount()];
                Point other_3 = new_other.Vertices[(i + ind + 2) % new_other.verticesCount()];

                Point other_reversed_1 = new_other.Vertices[save - (i + ind) % new_other.verticesCount()];
                Point other_reversed_2 = new_other.Vertices[save - (i + ind + 1) % new_other.verticesCount()];
                Point other_reversed_3 = new_other.Vertices[save - (i + ind + 2) % new_other.verticesCount()];

                Vector a = Vector(this_1, this_2);
                Vector a2 = Vector(this_2, this_3);
                Vector b, b2;
                if (!reversed_order) {
                    b = Vector(other_1, other_2);
                    b2 = Vector(other_2, other_3);
                } else {
                    b = Vector(other_reversed_1, other_reversed_2);
                    b2 = Vector(other_reversed_2, other_reversed_3);
                }
                if (std::abs(a.size() - b.size()) > EPS || std::abs(vector_product(a, a2) - vector_product(b, b2)) > EPS) {
                    fl = false;
                    break;
                }
            }
            if (fl)
                return true;
        }
        return false;
    }

public:
    bool isCongruentTo(const Shape& another) const override {
        try {
            const auto& other = dynamic_cast<const Polygon&>(another);

            std::vector<Point> new1 = delete_on_one_line(Vertices);
            std::vector<Point> new2 = delete_on_one_line(other.Vertices);
            Polygon new_this = Polygon(new1);
            Polygon new_other = Polygon(new2);
            if (new_this.verticesCount() != new_other.verticesCount()) {
                return false;
            }

            if (check_congruency(new_this, new_other, false))
                return true;
            if (check_congruency(new_this, new_other, true))
                return true;

            new_other.reflex(Line(0, 0));

            if (check_congruency(new_this, new_other, false))
                return true;
            if (check_congruency(new_this, new_other, true))
                return true;

            return false;
        } catch (std::bad_cast& b) {
            return false;
        }
        return false;
    }

    bool isSimilarTo(const Shape& another) const override {
        try {
            const auto& other = dynamic_cast<const Polygon&>(another);

            std::vector<Point> new1 = delete_on_one_line(Vertices);
            std::vector<Point> new2 = delete_on_one_line(other.Vertices);
            Polygon new_this = Polygon(new1);
            Polygon new_other = Polygon(new2);
            if (new_this.verticesCount() != new_other.verticesCount())
                return false;

            Vector v(new_this.Vertices[0], new_this.Vertices[1]);
            Vector v2(new_other.Vertices[0], new_other.Vertices[1]);
            long double mini1 = v.size();
            long double mini2 = v2.size();
            for (int i = 0; i < new_this.verticesCount(); ++i) {
                v = Vector(new_this.Vertices[i], new_this.Vertices[(i + 1) % new_this.verticesCount()]);
                mini1 = std::min(mini1, v.size());
            }
            for (int i = 0; i < new_other.verticesCount(); ++i) {
                v2 = Vector(new_other.Vertices[i], new_other.Vertices[(i + 1) % new_other.verticesCount()]);
                mini2 = std::min(mini2, v2.size());
            }
            new_other.scale(Point(0, 0), mini1 / mini2);
            return new_this.isCongruentTo(new_other);
        } catch (std::bad_cast& b) {
            return false;
        }
        return false;
    }

    bool containsPoint(const Point& point) const override {
        long double k = 7;
        Line check(point, k);
        int cnt = 0;
        for (int i = 0; i < verticesCount(); ++i) {
            Line l(Vertices[i], Vertices[(i + 1) % verticesCount()]);
            Point m = l.intersection_point(check);
            Vector fi(m, Vertices[i]);
            Vector se(m, Vertices[(i + 1) % verticesCount()]);
            if (m.x >= point.x - EPS && scalar_product(fi, se) <= 0) {
                cnt++;
            }
        }
        return cnt % 2 == 1;
    }

    void rotate(const Point& center, long double angle) override {
        for (auto& point : Vertices) {
            Vector a = Vector(center, point);
            a.rotate(angle);
            point = center + a;
        }
    }

    void reflex(const Point& center) override {
        for (auto& point : Vertices) {
            Vector to_center = Vector(point, center);
            point = center + to_center;
        }
    }

    void reflex(const Line& axis) override {
        for (auto& point : Vertices) {
            long double hight = axis.hight(point);
            Vector n = axis.normal_vector();
            n /= n.size();
            n *= hight;
            if (std::abs(axis.result(point + n)) > EPS)
                n *= -1;
            n *= 2;
            point = point + n;
        }
    }

    void scale(const Point& center, long double coefficient) override {
        for (auto& point : Vertices) {
            Vector from_center = Vector(center, point);
            from_center *= coefficient;
            point = center + from_center;
        }
    }
};

class Ellipse: public Shape {
protected:
    Point Focus1;
    Point Focus2;
    long double A;
    long double C;
    long double B;

public:
    Ellipse(Point a, Point b, long double sum): Focus1(a), Focus2(b),
        A(sum / 2), C(space(a, b) / 2), B(sqrtl(A * A - C * C)) {}

    std::pair<Point,Point> focuses() const {
        std::pair<Point,Point> ans = {Focus1, Focus2};
        return ans;
    }

    std::pair<Line, Line> directrices() const {
        Vector v(Focus1, Focus2);
        v /= v.size();
        v *= (A * A / C);
        Point p1 = center() + v;
        v *= -1;
        Point p2 = center() + v;
        Vector new_v(-v.y, v.x);
        std::pair<Line, Line> l = {Line(p1, p1 + new_v), Line(p2, p2 + new_v)};
        return l;
    }

    long double eccentricity() const {
        return C / A;
    }

    Point center() const {
        Vector a = Vector(Focus1, Focus2);
        a /= 2;
        Point ans = Focus1 + a;
        return ans;
    }

    long double perimeter() const override  {
        long double ans = PI * (3 * (A + B) - sqrtl((3 * A + B) * (A + 3 * B)));
        return ans;
    }

    long double area() const override {
        return PI * A * B;
    }

    bool operator== (const Shape& another) const override {
        try {
            const auto& other = dynamic_cast<const Ellipse&>(another);
            return std::abs(A - other.A) < EPS && std::abs(B - other.B) < EPS && std::abs(C - other.C) < EPS &&
                   ((Focus1 == other.Focus1 && Focus2 == other.Focus2) ||
                    (Focus1 == other.Focus2 && Focus2 == other.Focus1));
        } catch (std::bad_cast& b) {
            return false;
        }
        return false;
    }

    bool operator!= (const Shape& another) const override {
        return !(*this == another);
    }

    bool isCongruentTo(const Shape& another) const override {
        try {
            const auto& other = dynamic_cast<const Ellipse&>(another);
            return (std::abs(A - other.A) < EPS && std::abs(B - other.B) < EPS && std::abs(C - other.C) < EPS &&
                    std::abs(space(Focus1, Focus2) - space(other.Focus1, other.Focus2)) < EPS);
        } catch (std::bad_cast& b) {
            return false;
        }
        return false;
    }

    bool isSimilarTo(const Shape& another) const override {
        try {
            const auto& other = dynamic_cast<const Ellipse&>(another);
            long double coef = A / other.A;
            return (std::abs(B / other.B - coef) < EPS && std::abs(C / other.C - coef) < EPS &&
                    std::abs(space(Focus1, Focus2) / space(other.Focus1, other.Focus2) - coef) < EPS);
        } catch (std::bad_cast& b) {
            return false;
        }
        return false;
    }

    bool containsPoint(const Point& point) const override {
        return space(Focus1, point) + space(Focus2, point) < 2 * A + EPS;
    }

    void rotate(const Point& center, long double angle) override {
        Vector v(center, Focus1);
        v.rotate(angle);
        Focus1 = center + v;
        Vector v2(center, Focus2);
        v2.rotate(angle);
        Focus2 = center + v2;
    }

    void reflex(const Point& center) override {
        Vector to_center1 = Vector(Focus1, center);
        Focus1 = center + to_center1;
        Vector to_center2 = Vector(Focus2, center);
        Focus2 = center + to_center2;
    }

    void reflex(const Line& axis) override {
        long double hight = axis.hight(Focus1);
        Vector n = axis.normal_vector();
        n /= n.size();
        n *= hight;
        if (std::abs(axis.result(Focus1 + n)) > EPS)
            n *= -1;
        n *= 2;
        Focus1 = Focus1 + n;

        hight = axis.hight(Focus2);
        n = axis.normal_vector();
        n /= n.size();
        n *= hight;
        if (std::abs(axis.result(Focus2 + n)) > EPS)
            n *= -1;
        n *= 2;
        Focus2 = Focus2 + n;
    }

    void scale(const Point& center, long double coefficient) override {
        Vector v = Vector(center, Focus1);
        v *= coefficient;
        Focus1 = center + v;
        v = Vector(center, Focus2);
        v *= coefficient;
        Focus2 = center + v;
        A *= coefficient;
        B *= coefficient;
        C *= coefficient;
    }
};

class Circle: public Ellipse {
public:
    Circle(Point a, long double r): Ellipse(a, a, 2 * r) {}

    long double radius() const {
        return A;
    }

    Point center() const {
        return Focus1;
    }
};

class Rectangle: public Polygon {
public:
    Rectangle(std::vector<Point>& object): Polygon(object) {}

    Rectangle(std::initializer_list<Point> list): Polygon(list) {}

    Rectangle(Point a, Point b, long double c): Polygon({a, a, b, b}) {
        Vector v = Vector(a, b);
        v *= (c * c) / (c * c + 1);
        Point h = a + v;
        Vector v2 = Vector(-v.y, v.x);
        Point ans = h + v2;
        Vector ac = Vector(a, ans);
        ac *= -1;
        Point ans2 = b + ac;
        Vertices[1] = ans;
        Vertices[3] = ans2;
    }

    Point center() const {
        Vector v(Vertices[0], Vertices[2]);
        v /= 2;
        return Vertices[0] + v;
    }

    std::pair<Line, Line> diagonals() const {
        std::pair<Line, Line> l = {Line(Vertices[0], Vertices[2]), Line(Vertices[1], Vertices[3])};
        return l;
    }
};

class Square: public Rectangle {
public:
    Square(std::vector<Point>& object): Rectangle(object) {}

    Square(std::initializer_list<Point> list): Rectangle(list) {}

    Square(Point a, Point b): Rectangle(a, b, 1) {}

    Circle circumscribedCircle() const {
        Vector v(Vertices[0], Vertices[2]);
        Circle ans(center(), v.size() / 2);
        return ans;
    }

    Circle inscribedCircle() const {
        Vector v(Vertices[0], Vertices[1]);
        Circle ans(center(), v.size() / 2);
        return ans;
    }

};

class Triangle: public Polygon {
private:
    Point a;
    Point b;
    Point c;
public:
    Triangle(std::vector<Point>& object): Polygon(object), a(object[0]), b(object[1]), c(object[2]) {}

    Triangle(std::initializer_list<Point> list): Polygon(list) {
        a = Vertices[0];
        b = Vertices[1];
        c = Vertices[2];
    }

    Triangle(const Point& p0, const Point& p1, const Point& p2): Triangle({p0, p1, p2}) {}

    Circle circumscribedCircle() const {
        long double x12, x23, x31, y12, y23, y31;
        long double z1, z2, z3, zx, zy, z;
        x12 = a.x - b.x;
        x23 = b.x - c.x;
        x31 = c.x - a.x;
        y12 = a.y - b.y;
        y23 = b.y - c.y;
        y31 = c.y - a.y;
        z1 = (a.x * a.x) + (a.y * a.y);
        z2 = (b.x * b.x) + (b.y * b.y);
        z3 = (c.x * c.x) + (c.y * c.y);
        zx = y12 * z3 + y23 * z1 + y31 * z2;
        zy = x12 * z3 + x23 * z1 + x31 * z2;
        z = x12 * y31 - y12 * x31;
        Point ans(-zx / 2 / z, zy / 2 / z);
        Circle res(ans, space(ans, a));
        return res;
    }

    Circle inscribedCircle() const {
        long double as, bs, cs;
        as = Vector(b, c).size();
        bs = Vector(a, c).size();
        cs = Vector(a, b).size();
        Point ans = Point((as * a.x + bs * b.x + cs * c.x) / (as + bs + cs), (as * a.y + bs * b.y + cs * c.y) / (as + bs + cs));
        Line l(a, b);
        Circle res(ans, l.hight(ans));
        return res;
    }

    Point centroid() const {
        Point ans((a.x + b.x + c.x) / 3, (a.y + b.y + c.y) / 3);
        return ans;
    }

    Point orthocenter() const {
        Vector v(a, b);
        Vector v2(-v.y, v.x);
        Line fi(c, c + v2);

        v = Vector(a, c);
        v2 = Vector(-v.y, v.x);
        Line se(b, b + v2);

        Point ans = fi.intersection_point(se);
        return ans;
    }

    Line EulerLine() const {
        Circle res = circumscribedCircle();
        Line ans(res.center(), orthocenter());
        return ans;
    }

    Circle ninePointsCircle() const {
        Point new_a = middle(b, c);
        Point new_b = middle(a, c);
        Point new_c = middle(a, b);
        Triangle new_s({new_a, new_b, new_c});
        return new_s.circumscribedCircle();
    }
};