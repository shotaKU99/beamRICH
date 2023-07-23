#include "Vector3D.hh"

Vector3D::Vector3D() : x(1.0), y(0.0), z(0.0) {}

Vector3D::Vector3D(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}

Vector3D Vector3D::operator=(const Vector3D& right) { 
    x = right.x;
    y = right.y;
    z = right.z;
    return (*this); 
}

Vector3D Vector3D::operator+(const Vector3D& right) {
    Vector3D v;
    v.x = x + right.x;
    v.y = y + right.y;
    v.z = z + right.z;
    return v;
}

Vector3D Vector3D::operator-(const Vector3D& right) {
    Vector3D v;
    v.x = x - right.x;
    v.y = y - right.y;
    v.z = z - right.z;
    return v;
}

Vector3D Vector3D::operator-() {
    Vector3D v;
    v.x = -x;
    v.y = -y;
    v.z = -z;
    return v;
}

Vector3D Vector3D::operator*(const Vector3D& right) {
    Vector3D v;
    v.x = x * right.x;
    v.y = y * right.y;
    v.z = z * right.z;
    return v;
}

Vector3D Vector3D::operator/(const Vector3D& right) {
    Vector3D v;
    v.x = x / right.x;
    v.y = y / right.y;
    v.z = z / right.z;
    return v;
}

Vector3D& Vector3D::operator+=(const Vector3D& right) {
    x += right.x;
    y += right.y;
    z += right.z;
    return (*this);
}

Vector3D& Vector3D::operator-=(const Vector3D& right) {
    x -= right.x;
    y -= right.y;
    z -= right.z;
    return (*this);
}
Vector3D& Vector3D::operator*=(const Vector3D& right) {
    x *= right.x;
    y *= right.y;
    z *= right.z;
    return (*this);
}
Vector3D& Vector3D::operator/=(const Vector3D& right) {
    x /= right.x;
    y /= right.y;
    z /= right.z;
    return (*this);
}

Vector3D& Vector3D::operator+=(double a) {
    x += a;
    y += a;
    z += a;
    return (*this);
}
Vector3D& Vector3D::operator-=(double a) {
    x -= a;
    y -= a;
    z -= a;
    return (*this);
}
Vector3D& Vector3D::operator*=(double a) {
    x *= a;
    y *= a;
    z *= a;
    return (*this);
}
Vector3D& Vector3D::operator/=(double a) {
    x /= a;
    y /= a;
    z /= a;
    return (*this);
}

double Vector3D::norm() const { return std::sqrt(norm2()); }

double Vector3D::norm2() const { return x * x + y * y + z * z; }

void Vector3D::normalize() { *this = normalized(); }

Vector3D Vector3D::normalized() const {
    double inv_norm = 1.0 / norm();
    return Vector3D(x, y, z) * inv_norm;
}

Vector3D Vector3D::rotateX(double angleX) {
    return Vector3D(x, y * std::cos(angleX) - z * std::sin(angleX),
                    y * std::sin(angleX) + z * std::cos(angleX));
}

Vector3D Vector3D::rotateY(double angleY) {
    return Vector3D(x * std::cos(angleY) + z * std::sin(angleY), y,
                    - x * std::sin(angleY) + z * std::cos(angleY));
}

Vector3D Vector3D::rotateZ(double angleZ) {
    return Vector3D(x * std::cos(angleZ) - y * std::sin(angleZ),
                    x * std::sin(angleZ) + y * std::cos(angleZ), z);
}

bool Vector3D::isZero() { return (x == 0 && y == 0 && z == 0); }

double Vector3D::dot(const Vector3D left, const Vector3D right) {
    return left.x * right.x + left.y * right.y + left.z * right.z;
}

Vector3D Vector3D::cross(const Vector3D left, const Vector3D right) {
    return Vector3D(left.y * right.z - left.z * right.y, left.z * right.x - left.x * right.z,
                    left.x * right.y - left.y * right.x);
}

Vector3D Vector3D::reflect(const Vector3D incident, const Vector3D nperp) {
    Vector3D inci_norm = incident.normalized();
    Vector3D nperp_norm = nperp.normalized();

    return (inci_norm - 2.0 * dot(inci_norm, nperp_norm) * nperp_norm);
}

Vector3D Vector3D::refraction(const Vector3D incident, const Vector3D nperp, const double n1,
                              const double n2) {
    Vector3D inci_norm = incident.normalized();
    Vector3D nperp_norm = nperp.normalized();
    double costheta = std::abs(dot(inci_norm, nperp_norm));
    double k = 1.0 - std::pow(n1 / n2, 2.0) * (1.0 - std::pow(costheta, 2.0));
    if (k < 0.0) {
        //std::cout << "Refraction error occured" << std::endl;
        //std::cout << "costheta = " << costheta << std::endl;
        //std::cout << "k = " << k << std::endl;
        return Vector3D(0.0, 0.0, 1.0);
    } else {
        return n1 * inci_norm / n2 + (std::sqrt(k) - n1 * costheta / n2) * nperp_norm;
    }
}

std::ostream& operator<<(std::ostream& os, const Vector3D& vector3d) {
    os << '{' << vector3d.x << ',' << vector3d.y << ',' << vector3d.z << '}';
    return os;
}

std::istream& operator>>(std::istream& is, Vector3D& vector3d) {
    is.ignore(256, '{');
    is >> vector3d.x;
    is.ignore(256, ',');
    is >> vector3d.y;
    is.ignore(256, ',');
    is >> vector3d.z;
    is.ignore(256, '}');
    return is;
}

Vector3D operator+(double a, const Vector3D& right) {
    Vector3D v;
    v.x = right.x + a;
    v.y = right.y + a;
    v.z = right.z + a;
    return v;
};
Vector3D operator+(const Vector3D& left, double a) {
    Vector3D v;
    v.x = left.x + a;
    v.y = left.y + a;
    v.z = left.z + a;
    return v;
};

Vector3D operator-(double a, const Vector3D& right) {
    Vector3D v;
    v.x = right.x - a;
    v.y = right.y - a;
    v.z = right.z - a;
    return v;
};
Vector3D operator-(const Vector3D& left, double a) {
    Vector3D v;
    v.x = left.x - a;
    v.y = left.y - a;
    v.z = left.z - a;
    return v;
};

Vector3D operator*(double a, const Vector3D& right) {
    Vector3D v;
    v.x = right.x * a;
    v.y = right.y * a;
    v.z = right.z * a;
    return v;
};
Vector3D operator*(const Vector3D& left, double a) {
    Vector3D v;
    v.x = left.x * a;
    v.y = left.y * a;
    v.z = left.z * a;
    return v;
};
Vector3D operator/(const Vector3D& left, double a) {
    Vector3D v;
    v.x = left.x / a;
    v.y = left.y / a;
    v.z = left.z / a;
    return v;
};
