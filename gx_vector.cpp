#include "gx_vector.h"
#include <iostream>

static inline double min3(double a, double b, double c)
{
    return std::min(a, std::min(b,c));
}

static inline double max3(double a, double b, double c)
{
    return std::max(a, std::max(b,c));
}

/* Constructors */
Vector::Vector(){
    coords[0] = 0.;
    coords[1] = 0.;
    coords[2] = 0.;
}
Vector::Vector(double x,double y, double z){
    coords[0] = x;
    coords[1] = y;
    coords[2] = z;
}

Vector::Vector(const Vector& other)
{   
    for(int i=0; i<3; i++)
        coords[i] = other[i];
}

Vector::Vector(Vector&& other)
{
    for(int i=0; i<3; i++)
        coords[i] = other[i];
}

/* Getters and Setters*/
double Vector::operator[](int i) const {
    return coords[i];
}

double& Vector::operator[](int i){
    return coords[i];
}

/* Inplace Operators */
Vector& Vector::operator+=(const Vector& rhs){
    for(int i=0; i<3; i++)
        coords[i] += rhs[i];
    return *this;
}
Vector& Vector::operator-=(const Vector& rhs){
    for(int i=0; i<3; i++)
        coords[i] -= rhs[i];
    return *this;
}
Vector& Vector::operator*=(const Vector& rhs){
    for(int i=0; i<3; i++)
        coords[i] *= rhs[i];
    return *this;
}
Vector& Vector::operator*=(double rhs){
    for(int i=0; i<3; i++)
        coords[i] *= rhs;
    return *this;
}
Vector& Vector::operator/=(const Vector& rhs){
    for(int i=0; i<3; i++)
        coords[i] /= rhs[i];
    return *this;
}
Vector& Vector::operator/=(double rhs){
    for(int i=0; i<3; i++)
        coords[i] /= rhs;
    return *this;
}
Vector& Vector::operator=(const Vector& rhs){
    for(int i=0; i<3; i++)
        coords[i] = rhs[i];
    return *this;
}

// Not inplace operators 
Vector operator+(const Vector& a, const Vector& b){
    return Vector(a[0]+b[0], a[1]+b[1], a[2]+b[2]);
}

Vector operator-(const Vector& a, const Vector& b) {
	return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const Vector& a, const Vector& b){
	/* coordinate-wise product */
	return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
Vector operator*(const double a, const Vector& b) {
	return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
	return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator/(const Vector& a, const Vector& b){
    return Vector(a[0]/b[0], a[1]/b[1], a[2]/b[2]);
}
Vector operator/(const Vector& a, const double b) {
	return (1/b) * a;
}
Vector operator-(const Vector& a) {
	return Vector(-a[0], -a[1], -a[2]);
}

double dot(const Vector& a, const Vector& b){
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
double cosine_similarity(const Vector& a, const Vector& b){
	double prod_norm = a.norm() * b.norm();
	double dot_prod = dot(a, b);
	return dot_prod/prod_norm;
}
Vector cross(const Vector& a, const Vector& b) {
	return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}

/* Utils */

double Vector::min() const
{
    return min3(coords[0], coords[1], coords[2]);
}

double Vector::max() const
{
    return max3(coords[0], coords[1], coords[2]);
}

// Norm
double Vector::norm2() const{
    return (
        coords[0] * coords[0] + 
        coords[1] * coords[1] + 
        coords[2] * coords[2]
    );
}
double Vector::norm() const{
    return sqrt(norm2());
}
Vector Vector::normalized() const{
    return *this/norm();
}

// Other operators
Vector Vector::abs() const{
    return Vector(
        std::abs(coords[0]),
        std::abs(coords[1]),
        std::abs(coords[2])
    );
}
Vector Vector::clip(double min_value, double max_value) const{
    Vector output(*this);
    for(int i=0; i<3; i++)
        output[i] = std::max(std::min(output[i], max_value), min_value);
    return output;
}

// Coordinates
size_t Vector::argmin() const{
    if(coords[0] <= coords[1] && coords[0] <= coords[2])
        return 0;
    else if(coords[1] <= coords[2])
        return 1;
    return 2;
}
size_t Vector::argmax() const{
    if(coords[0] >= coords[1] && coords[0] >= coords[2])
        return 0;
    else if(coords[1] >= coords[2])
        return 1;
    return 2;
}

// Other
Vector Vector::orthogonal() const{
    Vector x = *this;
    size_t i = x.abs().argmin();
    x[i] = 0.;
	if(i==0){
		x[1] = -x[1];
		std::swap(x[1], x[2]);
	} 
	else if(i==1){
		x[2] = -x[2];
		std::swap(x[2], x[0]);
	}
	else{
		x[0] = -x[0];
		std::swap(x[0], x[1]);
	}
	return x;
}

/* Print */
std::ostream& operator<<(std::ostream& os, const Vector& rhs){
    os << "Vector(" << rhs[0] << "," << rhs[1] << "," << rhs[2] << ")";
    return os;
}

// static Vector Methods

Vector Vector::min(const Vector& a, const Vector& b)
{
    return Vector(std::min(a[0], b[0]), std::min(a[1], b[1]), std::min(a[2], b[2]));
}

Vector Vector::max(const Vector& a, const Vector& b)
{
    return Vector(std::max(a[0], b[0]), std::max(a[1], b[1]), std::max(a[2], b[2]));
}

Vector Vector::initOnesVector()
{
    return Vector(1.,1.,1.);
}

Vector Vector::initMinVector()
{
    return Vector(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());
}

Vector Vector::initMaxVector()
{
    return Vector(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
}

// ----- Ray Implementation -----
Ray::Ray(const Vector& O, const Vector& u)
: O(O), u(u)
{}

// ----- Matrix Implementation
Matrix::Matrix(const Vector& q1, const Vector& q2, const Vector& q3)
{
    q[0] = q1;
    q[1] = q2;
    q[2] = q3;
}

Matrix::Matrix(const double matrix[3][3])
{
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            q[i][j] = matrix[i][j];
}

Vector Matrix::operator[](int i) const
{
    return q[i];
}

Vector& Matrix::operator[](int i)
{
    return q[i];
}

Matrix& Matrix::operator*=(double rhs)
{
    for(auto &vec : q)
        vec*=rhs;
    return *this;
}

Matrix& Matrix::operator/=(double rhs)
{
    return (*this) *= 1/rhs;
}

Matrix& Matrix::operator*=(const Matrix& rhs)
{
    double out[3][3] = {};
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            out[i][j] = 0.;

    for(int k=0; k<3; k++)
        for(int j=0; j<3; j++)
            for(int i=0; i<3; i++)
                out[i][j] += q[i][k] * rhs[k][j];
    
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            q[i][j] = out[i][j];
    
    return *this;
}

// Not inplace operators
Matrix operator+(const Matrix& a, const Matrix& b)
{
    return Matrix(a[0]+b[0], a[1]+b[1], a[2]+b[2]);
}
Matrix operator-(const Matrix& a, const Matrix& b)
{
    return Matrix(a[0]-b[0], a[1]-b[1], a[2]-b[2]);
}
Matrix operator*(const Matrix& a, const Matrix& b)
{
    double out[3][3] = {};
    for(int i=0; i<3; i++)
        for(int j=0; j<3; j++)
            out[i][j] = 0.;

    for(int k=0; k<3; k++)
        for(int j=0; j<3; j++)
            for(int i=0; i<3; i++)
                out[i][j] += a[i][k] * b[k][j];
    
    return Matrix(out);
}

Matrix operator*(const double a, const Matrix& b)
{
    return Matrix(a*b[0], a*b[1], a*b[2]);
}
Matrix operator*(const Matrix& a, const double b)
{
    return Matrix(a[0]*b, a[1]*b, a[2]*b);
}
Matrix operator/(const Matrix& a, const double b)
{
    return a*(1/b);
}
Matrix operator-(const Matrix& a)
{
    return Matrix(-a[0], -a[1], -a[2]);
}
Vector operator*(const Matrix& a, const Vector& b)
{
    return Vector(
        dot(a[0],b),
        dot(a[1],b),
        dot(a[2],b)
    );
}

std::ostream& operator<<(std::ostream& os, const Matrix& rhs)
{
    os << "Matrix("<<rhs[0]<<", "<<rhs[1]<<", "<<rhs[2]<<")";
    return os;
}