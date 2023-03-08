#ifndef VECTOR3D_h
#define VECTOR3D_h

#include <iostream>
#include <cmath>

class Vector3D{

public:
    double x, y, z;
    Vector3D();
    Vector3D(const double _x, const double _y, const double _z);

    double operator[](int index) const{
        if(index==0){
            return x;
        }else if(index == 1){
            return y;
        }else if(index == 2){
            return z;
        }else{
            printf("Index Error!!");
            exit (1);
        }
    };

    // Vector - Vector operators
    Vector3D operator+(const Vector3D &right);
    Vector3D operator-(const Vector3D &right);
    Vector3D operator-();
    Vector3D operator*(const Vector3D &right);
    Vector3D operator/(const Vector3D &right);
    Vector3D& operator+=(const Vector3D &right);
    Vector3D& operator-=(const Vector3D &right);
    Vector3D& operator*=(const Vector3D &right);
    Vector3D& operator/=(const Vector3D &right);

    // Vector - Scalar operators
    Vector3D& operator+=(const double a);
    Vector3D& operator-=(const double a);
    Vector3D& operator*=(const double a);
    Vector3D& operator/=(const double a);

    double norm() const;
    double norm2() const;
    void normalize();
    Vector3D normalized() const;
    bool isZero();
    static double dot(const Vector3D& left, const Vector3D& right);
    static Vector3D cross(const Vector3D& left, const Vector3D& right);

    static Vector3D reflect(const Vector3D& incident, const Vector3D& nperp);
    static Vector3D refraction(const Vector3D& incident, const Vector3D& nperp, const double n1, const double n2);
    
};

std::ostream& operator<<(std::ostream& os, const Vector3D& vector3);
std::istream& operator>>(std::istream& is, Vector3D& vector3);

// Vector - Scalar operators
Vector3D operator+(const double a, const Vector3D& right);
Vector3D operator+(const Vector3D& left, const double a);

Vector3D operator-(const double a, const Vector3D& right);
Vector3D operator-(const Vector3D& left, const double a);

Vector3D operator*(const double a, const Vector3D& right);
Vector3D operator*(const Vector3D& left, const double a);

Vector3D operator/(const Vector3D& left, const double a);




#endif