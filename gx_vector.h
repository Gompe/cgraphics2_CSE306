#pragma once

#include <ostream>
#include <vector>
#include <algorithm>
#include <cmath>

// Vector Declaration
class Vector{
public:
    Vector();
    explicit Vector(double x,double y, double z);
    Vector(const Vector& other);
    Vector(Vector&& other);

    double operator[](int i) const;
    double& operator[](int i);

    // Inplace operators
    Vector& operator+=(const Vector& rhs);
    Vector& operator-=(const Vector& rhs);
    Vector& operator*=(const Vector& rhs);
    Vector& operator*=(double rhs);
    Vector& operator/=(const Vector& rhs);
    Vector& operator/=(double rhs);
    Vector& operator=(const Vector& rhs);

    double min() const;
    double max() const;

    double norm2() const;
    double norm() const;
    Vector normalized() const;
    Vector abs() const;
    Vector clip(double min_value, double max_value) const;
    size_t argmin() const;
    size_t argmax() const;
    Vector orthogonal() const;

    static Vector min(const Vector& a, const Vector& b);
    static Vector max(const Vector& a, const Vector& b);

    static Vector initOnesVector();
    static Vector initMinVector();
    static Vector initMaxVector();


private:
    double coords[3];
};

// Not inplace operators
Vector operator+(const Vector& a, const Vector& b);
Vector operator-(const Vector& a, const Vector& b);
Vector operator*(const Vector& a, const Vector& b);
Vector operator*(const double a, const Vector& b);
Vector operator*(const Vector& a, const double b);
Vector operator/(const Vector& a, const Vector& b);
Vector operator/(const Vector& a, const double b);
Vector operator-(const Vector& a);

double dot(const Vector& a, const Vector& b);
double cosine_similarity(const Vector& a, const Vector& b);
Vector cross(const Vector& a, const Vector& b);

std::ostream& operator<<(std::ostream& os, const Vector& rhs);

// Ray Declaration
class Ray
{
public:
    Ray(const Vector&O, const Vector& u);
    Vector O;
    Vector u;
};

// Declaration of 3x3 Matrix
class Matrix
{
public:
    Matrix() = default;

    // Creates matrix row by row
    explicit Matrix(const Vector& q1, const Vector& q2, const Vector& q3);
    explicit Matrix(const double matrix[3][3]);

    // Gets row i
    Vector operator[](int i) const;
    Vector& operator[](int i);

    Matrix& operator*=(double rhs);
    Matrix& operator/=(double rhs);
    Matrix& operator*=(const Matrix& rhs);

    Matrix transpose() const;
    Matrix& inplaceTranspose();

private:
    Vector q[3];
};

Matrix operator+(const Matrix& a, const Matrix& b);
Matrix operator-(const Matrix& a, const Matrix& b);
Matrix operator*(const Matrix& a, const Matrix& b);
Matrix operator*(const double a, const Matrix& b);
Matrix operator*(const Matrix& a, const double b);
Matrix operator/(const Matrix& a, const double b);
Matrix operator-(const Matrix& a);
Vector operator*(const Matrix& a, const Vector& b);

std::ostream& operator<<(std::ostream& os, const Matrix& rhs);
