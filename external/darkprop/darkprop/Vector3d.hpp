/**
 * @file
 * We just need a simple Euclidean 3 vector.
 */
#pragma once
#include <cmath>
#include <string>
#include <vector>
#include <ostream>
#include <stdexcept>

namespace darkprop {

template<typename Value>
class Vector3d
{
public:
    Value x, y, z;
    Vector3d() {}
    Vector3d(Value t_x, Value t_y, Value t_z) : x {t_x}, y {t_y}, z {t_z} {}
    explicit Vector3d(const std::vector<Value>& v) {
        if (v.size() != 3) { throw std::invalid_argument("length must be 3"); }
        x = v[0]; y = v[1]; z = v[2];
    }
    Vector3d(const Vector3d<Value>& other) :
        x {other.x}, y {other.y}, z {other.z} {}

    Vector3d<Value>& operator=(const Vector3d<Value>& other) {
        if (this == &other) return *this;
        x = other.x;
        y = other.y;
        z = other.z;
        return *this;
    }
    Vector3d<Value> operator+(const Vector3d<Value>& other) const {
        return Vector3d(x + other.x, y + other.y, z + other.z);
    }
    Vector3d<Value> operator+(Value s) const {
        return Vector3d(x + s, y + s, z + s);
    }
    Vector3d<Value> operator-(const Vector3d<Value>& other) const {
        return Vector3d(x - other.x, y - other.y, z - other.z);
    }
    Vector3d<Value> operator-(Value s) const {
        return Vector3d(x - s, y - s, z - s);
    }
    Vector3d<Value> operator-() const {
        return Vector3d(-x, -y, -z);
    }
    Vector3d<Value> operator/(Value s) const {
        return Vector3d(x / s, y / s, z / s);
    }
    Vector3d<Value> operator*(Value s) const {
        return Vector3d(x * s, y * s, z * s);
    }
    Vector3d<Value>& operator+=(const Vector3d<Value>& rhs) {
        x += rhs.x; y += rhs.y; z += rhs.z;
        return *this;
    }
    Vector3d<Value>& operator+=(Value s) {
        x += s; y += s; z += s;
        return *this;
    }
    Vector3d<Value>& operator-=(const Vector3d<Value>& rhs) {
        x -= rhs.x; y -= rhs.y; z -= rhs.z;
        return *this;
    }
    Vector3d<Value>& operator-=(Value s) {
        x -= s; y -= s; z -= s;
        return *this;
    }
    bool operator==(const Vector3d<Value>& rhs) const {
        return x == rhs.x && y == rhs.y && z == rhs.z;
    }
    bool isApprox(const Vector3d<Value>& other, Value rel_prec) const {
        Value n = norm();
        Value othern = other.norm();
        return std::abs(n - othern) < std::max(n, othern) * rel_prec;
    }

    friend Vector3d<Value> operator*(Value s, const Vector3d<Value>& v) {
        return v * s;
    }
    friend Vector3d<Value> operator+(Value s, const Vector3d<Value>& v) {
        return Vector3d(s + v.x, s + v.y, s + v.z);
    }
    friend Vector3d<Value> operator-(Value s, const Vector3d<Value>& v) {
        return Vector3d(s - v.x, s - v.y, s - v.z);
    }

    friend std::ostream& operator<<(std::ostream& os, const Vector3d<Value>& v) {
        os << "[" << v.x << ", " << v.y << ", " << v.z << "]" << std::endl;
        return os;
    }

    Value norm() const { return std::hypot(x, y, z); }
    Vector3d<Value>& normalize() {
        Value n = norm();
        if (n != 0) {
            x = x / n; y = y / n; z = z / n;
        }
        return *this;
    }
    Vector3d<Value> normalized() const {
        return Vector3d<Value>(*this).normalize();
    }
    Value dot(const Vector3d<Value>& other) const {
        return x * other.x + y * other.y + z * other.z;
    }
    Vector3d<Value> cross(const Vector3d<Value>& other) const {
        return Vector3d(y * other.z - other.y * z,
                        z * other.x - other.z * x,
                        x * other.y - other.x * y);
    }

    const Value& operator[](int i) const {
        if (i == 0) return x;
        if (i == 1) return y;
        if (i == 2) return z;
        throw std::out_of_range {"index i = " + std::to_string(i)};
    }
    Value& operator[](int i) {
        return const_cast<Value&>(
                static_cast<const Vector3d<Value>&>(*this)[i]);
    }
    Value& operator[](std::size_t i) { return operator[](int(i)); }
    const Value& operator[](std::size_t i) const { return operator[](int(i)); }
    const Value& operator()(int i) const { return operator[](i); }
    const Value& operator()(std::size_t i) const { return operator[](i); }

    std::vector<Value> to_vec() const { return {x, y, z}; }
};

} // namespace darkprop
