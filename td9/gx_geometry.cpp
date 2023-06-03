#include "gx_geometry.h"

#include <iostream> 
#include<assert.h>
#include <functional>

template <typename T>
static inline void print(const T& x, bool endLine=true) {
    if (endLine)
        std::cout << x << std::endl;
    else
        std::cout << x << " ";
}

static inline double squareDouble(double x)
{
    return x*x;
}

static inline double max3(double a, double b, double c)
{
    return std::max(a, std::max(b, c));
}

static inline double min3(double a, double b, double c)
{
    return std::min(a, std::min(b, c));
}


Sphere::Sphere(const Vector& C, double R)
: C(C), R(R)
{
    initBbox();
}

Sphere::Sphere(const Vector& C, double R, const Material& m)
: Object(m), C(C), R(R)
{
    initBbox();
}

void Sphere::initBbox()
{   
    Vector ones(1., 1., 1.);
    bbox = BoundingBox(C-R*ones, C+R*ones);
}

bool Sphere::intersect(const Ray& ray, ObjectHit& hitInfo) const {
    Vector d = ray.O-C;
    double a = dot(ray.u, d);

    double delta = squareDouble(a)-d.norm2()+R*R;
    if (delta < 0)
        return false;
    double t2 = sqrt(delta) - a;
    if (t2 < 0)
        return false;
        
    double t1 = -sqrt(delta) - a;
    hitInfo.tHit = (t1 > 0) ? t1 : t2;
    hitInfo.P = ray.O + hitInfo.tHit*ray.u;
    hitInfo.N = (hitInfo.P - C).normalized();
    hitInfo.object_ptr = this;

    return true;
}

bool Sphere::updateIntersect(const Ray& ray, ObjectHit& hitInfo) const {
    Vector d = ray.O-C;
    double a = dot(ray.u, d);

    double delta = squareDouble(a)-d.norm2()+R*R;
    if (delta < 0)
        return false;

    double t2 = sqrt(delta) - a;
    if (t2 < 0)
        return false;
        
    double t1 = -sqrt(delta) - a;
    double t = (t1 > 0) ? t1 : t2;

    if (t >= hitInfo.tHit)
        return false;

    hitInfo.tHit = t;
    hitInfo.P = ray.O + hitInfo.tHit*ray.u;
    hitInfo.N = (hitInfo.P - C).normalized();
    hitInfo.object_ptr = this;

    return true;
}

// ----- Bounding Box ----- //
BoundingBox::BoundingBox(const Vector& lbBbox, const Vector& ubBbox)
: lbBbox(lbBbox), ubBbox(ubBbox)
{}

bool BoundingBox::intersect(const Ray& ray) const
{   
    Vector u = ray.u;
    const Vector& O = ray.O;

    // To avoid division by 0
    static const double eps = 1E-6;
    for (int i=0; i<3; i++) {
        if (std::abs(u[i]) < eps)
            u[i] = eps;
    }

    Vector t0, t1;
    t0 = (lbBbox - O)/u;
    t1 = (ubBbox - O)/u;

    for (int i=0; i<3; i++) {
        if (t1[i] < t0[i])
            std::swap(t0[i], t1[i]);
    }

    // Return true if t0 and t1 intersect on a positive value
    return (t1.min() >= 0) && (t1.min() >= t0.max());
}

bool BoundingBox::intersect(const Ray& ray, double& t_lower_bound) const
{   
    Vector u = ray.u;

    // To avoid division by 0
    const double eps = 1E-6;
    for (int i=0; i<3; i++) {
        if (std::abs(u[i]) < eps)
            u[i] = eps;
    }

    Vector t0 = (lbBbox - ray.O)/u;
    Vector t1 = (ubBbox - ray.O)/u;

    for (int i=0; i<3; i++) {
        if (t1[i] < t0[i])
            std::swap(t0[i], t1[i]);
    }

    t_lower_bound = std::max((double) 0, t0.max());

    // Return true if t0 and t1 intersect on a positive value
    return (t1.min() >= t_lower_bound);
}

bool BoundingBox::isInside(const Vector& v) const
{
    const double eps = 1E-6;

    return (
        (v[0] <= ubBbox[0] + eps) && 
        (v[1] <= ubBbox[1] + eps) &&
        (v[2] <= ubBbox[2] + eps) &&
        (v[0] >= lbBbox[0] - eps) && 
        (v[1] >= lbBbox[1] - eps) &&
        (v[2] >= lbBbox[2] - eps)
    );
}

void BoundingBox::resizeAddVector(const Vector& v)
{
    lbBbox = Vector::min(lbBbox, v);
    ubBbox = Vector::max(ubBbox, v);
}

void BoundingBox::merge(const BoundingBox& other)
{
    lbBbox = Vector::min(lbBbox, other.lbBbox);
    ubBbox = Vector::max(ubBbox, other.ubBbox);
}

BoundingBox BoundingBox::merge(const BoundingBox& a, const BoundingBox& b)
{
    return BoundingBox(
        Vector::min(a.lbBbox, b.lbBbox),
        Vector::max(a.ubBbox, b.ubBbox)
    );
}

BoundingBox BoundingBox::initWithLimits()
{   
    return BoundingBox(Vector::initMaxVector(), Vector::initMinVector());
}

// Transformations
static inline Matrix rotationMatrix(const Vector& axis, double theta)
{   
    // https://en.wikipedia.org/wiki/Rotation_matrix
    // Assumes axis is normalized
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);

    const Vector& u = axis;

    double out[3][3] = {
        {cos_theta+squareDouble(u[0])*(1-cos_theta), u[0]*u[1]*(1-cos_theta)-u[2]*sin_theta, u[0]*u[2]*(1-cos_theta)+u[1]*sin_theta},
        {u[1]*u[0]*(1-cos_theta)+u[2]*sin_theta, cos_theta+squareDouble(u[1])*(1-cos_theta), u[1]*u[2]*(1-cos_theta)-u[0]*sin_theta},
        {u[2]*u[0]*(1-cos_theta)-u[1]*sin_theta, u[2]*u[1]*(1-cos_theta)+u[0]*sin_theta, cos_theta+squareDouble(u[2])*(1-cos_theta)}
    };

    return Matrix(out);
}

static inline Matrix rotationMatrix(int axis, double theta)
{
    switch (axis) {
        case 1: return rotationMatrix(Vector(1,0,0), theta);
        case 2: return rotationMatrix(Vector(0,1,0), theta);
        case 3: return rotationMatrix(Vector(0,0,1), theta);
    }
    throw std::logic_error("Axis must be 0, 1, 2 in rotationMatrix");
}

void Sphere::transformTranslate(const Vector& delta)
{
    C+=delta;
}
void Sphere::transformScale(double r){
    R*=r;
    C*=r;
}
void Sphere::transformRotate(const Vector& axis, double theta)
{
    C=rotationMatrix(axis, theta)*C;
}
void Sphere::transformRotate(int axis, double theta)
{
    C=rotationMatrix(axis, theta)*C;
}

// ----- Triangle meshes ------
static bool intersectTriangle(const Ray& ray, const Vector& A, 
const Vector& B, const Vector& C, double data[4])
{   
    // Previous implementation
    Vector e1 = B-A;
    Vector e2 = C-A;
    Vector d = A - ray.O;
    const Vector& u = ray.u;
    
    Vector N = cross(e1, e2);
    double Z = 1/dot(u, N);

    double beta = dot(e2, cross(d, u)) * Z;
    double gamma = -dot(e1, cross(d, u)) * Z;
    double alpha = 1 - beta - gamma;
    double t = dot(d, N) * Z;

    if (alpha >= 0 && beta >= 0 && gamma >= 0 && t >= 0) {
        // Maybe change the way I do this with a struct
        data[0] = alpha;
        data[1] = beta;
        data[2] = gamma;
        data[3] = t;
        return true;
    }

    return false;
}

static inline BoundingBox triangleBoundingBox(const Vector& A, const Vector& B, const Vector& C)
{
    BoundingBox bbox = BoundingBox::initWithLimits();
    bbox.resizeAddVector(A);
    bbox.resizeAddVector(B);
    bbox.resizeAddVector(C);

    return bbox;
}


TriangleMesh::TriangleMesh(const char* obj)
{
    readOBJ(obj);

    if (vertices.size() == 0)
        throw std::logic_error("Triangle Mesh cannot be empty");

    buildBVH();

    bbox = BoundingBox::initWithLimits();
    for (const auto& v : vertices) {
        bbox.resizeAddVector(v);
    }
}

TriangleMesh::TriangleMesh(const char* obj, const Material& m)
: Object(m)
{
    readOBJ(obj);

    if (vertices.size() == 0)
        throw std::logic_error("Triangle Mesh cannot be empty");

    buildBVH();

    bbox = BoundingBox::initWithLimits();
    for (const auto& v : vertices) {
        bbox.resizeAddVector(v);
    }
}

bool TriangleMesh::intersect(const Ray& ray, ObjectHit& gHit) const
{
    throw std::logic_error("Not implemented error");
}

bool TriangleMesh::updateIntersect(const Ray& ray, ObjectHit& hitInfo) const
{   
    return fastUpdateIntersect(ray, hitInfo);

    double t_lower_bound;

    if (!bbox.intersect(ray, t_lower_bound)) {
        return false;
    }

    if (t_lower_bound > hitInfo.tHit)
        return false;

    bool hasUpdate = false;

    for(const auto& tIndex : indices) {
        double data[4];
        const Vector& A = vertices[tIndex.vtxi];
        const Vector& B = vertices[tIndex.vtxj];
        const Vector& C = vertices[tIndex.vtxk];
        if (intersectTriangle(ray, A, B, C, data)) {
            if (data[3] < hitInfo.tHit) {
                hasUpdate = true;
                hitInfo.tHit = data[3];
                hitInfo.P = data[0]*A + data[1]*B + data[2]*C;
                hitInfo.N = (data[0]*normals[tIndex.ni] + data[1]*normals[tIndex.nj] + data[2]*normals[tIndex.nk]).normalized();
                hitInfo.object_ptr = this;
            }
        }
    }

    return hasUpdate;
}

bool TriangleMesh::fastUpdateIntersect(const Ray& ray, ObjectHit& hitInfo) const
{
    return BVHIntersectRoutine(ray, hitInfo, treeBVH.root);
}

bool TriangleMesh::BVHIntersectRoutine(const Ray& ray, ObjectHit& hitInfo, const binary_tree::Node<dataBVH> *node) const
{
    assert(node != nullptr);

    double t_lower_bound;

    if (!node->data.bbox.intersect(ray, t_lower_bound)) 
        return false;

    if (t_lower_bound > hitInfo.tHit)
        return false;
    
    // Intersection is possible
    if (node->is_leaf()) {
        bool hasUpdate = false;

        for(int i=node->data.start; i < node->data.end; i++) {
            const auto& tIndex = indices[i];
            double data[4];
            const Vector& A = vertices[tIndex.vtxi];
            const Vector& B = vertices[tIndex.vtxj];
            const Vector& C = vertices[tIndex.vtxk];
            if (intersectTriangle(ray, A, B, C, data)) {
                if (data[3] < hitInfo.tHit) {
                    hasUpdate = true;
                    hitInfo.tHit = data[3];
                    hitInfo.P = data[0]*A + data[1]*B + data[2]*C;
                    hitInfo.N = (data[0]*normals[tIndex.ni] + data[1]*normals[tIndex.nj] + data[2]*normals[tIndex.nk]).normalized();
                    hitInfo.object_ptr = this;
                }
            }
        }
        return hasUpdate;
    } 
    else {
        bool resultLeft = BVHIntersectRoutine(ray, hitInfo, node->left);
        bool resultRight = BVHIntersectRoutine(ray, hitInfo, node->right);

        return resultLeft || resultRight;
    }
}

void TriangleMesh::transformBVH(std::function<void (BoundingBox&)> transform) 
{
    std::vector<binary_tree::Node<dataBVH> *> nodeStack;
    nodeStack.push_back(treeBVH.root);

    while (nodeStack.size()) {
        binary_tree::Node<dataBVH> *node = nodeStack.back();
        nodeStack.pop_back();
        
        if (node) {
            transform(node->data.bbox);
            nodeStack.push_back(node->left);
            nodeStack.push_back(node->right);
        }
    }
}

void TriangleMesh::transformTranslate(const Vector& delta) {
    for (auto& vertex : vertices)
        vertex+=delta;

    auto bboxTranslate = [delta](BoundingBox& bbox){
        bbox.lbBbox += delta;
        bbox.ubBbox += delta;
    };

    transformBVH(bboxTranslate);
}

void TriangleMesh::transformScale(double r){
    for (auto& vertex : vertices)
        vertex *= r;

    auto bboxScale = [r](BoundingBox& bbox){
        bbox.lbBbox *= r;
        bbox.ubBbox *= r;
    };

    transformBVH(bboxScale);
}

void TriangleMesh::transformRotate(const Vector& axis, double theta){
    Matrix Q = rotationMatrix(axis, theta);
    for (auto& vertex : vertices)
        vertex = Q*vertex;
    for (auto& normal : normals)
        normal = Q*normal;

    // More complicated, need to rebuild BVH
    buildBVH();
}

void TriangleMesh::transformRotate(int axis, double theta){
    Matrix Q = rotationMatrix(axis, theta);
    for (auto& vertex : vertices)
        vertex = Q*vertex;
    for (auto& normal : normals)
        normal = Q*normal;

    // More complicated, need to rebuild BVH
    buildBVH();
}

Vector TriangleMesh::triangleCentroid(int idx) const {
    const TriangleIndices& tIndex = indices[idx];
    return (vertices[tIndex.vtxi] + vertices[tIndex.vtxj] + vertices[tIndex.vtxk])/3;
}

// Acceleration - BVH
void TriangleMesh::_verifyTree(binary_tree::Node<dataBVH> *node) const
{
    if (!node) return;

    for (int ptr = node->data.start; ptr < node->data.end; ptr++) {
        const auto& tidx = indices[ptr];
        assert (node->data.bbox.isInside(vertices[tidx.vtxi]));
        assert (node->data.bbox.isInside(vertices[tidx.vtxj]));
        assert (node->data.bbox.isInside(vertices[tidx.vtxk]));
    }

    _verifyTree(node->left);
    _verifyTree(node->right);
}

void TriangleMesh::buildBVH() 
{   
    if (treeBVH.root) {
        treeBVH.reset();
    }

    binary_tree::Node<dataBVH> *new_root = buildBVH(0, indices.size());
    treeBVH.root = new_root;
}

binary_tree::Node<dataBVH>* TriangleMesh::buildBVH(int start, int end)
{
    // Compute bounding box
    BoundingBox bbox = BoundingBox::initWithLimits();
    for (int i=start; i<end; i++) {
        const auto& tidx = indices[i];
        bbox.merge(triangleBoundingBox(vertices[tidx.vtxi], vertices[tidx.vtxj], vertices[tidx.vtxk]));
    }

    // Create current Node
    dataBVH data;
    data.start = start;
    data.end = end;
    data.bbox = bbox;

    binary_tree::Node<dataBVH> *node = new binary_tree::Node<dataBVH>(data);

    // initial condition
    if (end - start <= BVH_STOP) {
        return node;
    }

    // Compute longest axis to perform partition
    const Vector& lbBbox = bbox.lbBbox;
    const Vector& ubBbox = bbox.ubBbox;

    int axis = (ubBbox - lbBbox).argmax();
    double mid = (ubBbox[axis] + lbBbox[axis])/2.;

    // Partition scheme
    int ptr = start;
    int pivot = start;

    for (ptr = start; ptr < end; ptr++) {
        if (triangleCentroid(ptr)[axis] < mid) {
            assert (start <= ptr && ptr < end && pivot >= start && pivot < end);
            std::swap(indices[ptr], indices[pivot]);
            pivot++;
        }
    }

    // pivot - start = Number of elements left to mid
    if (pivot == start || pivot == end)
        return node;

    // Will continue
    node->left = buildBVH(start, pivot);
    node->right = buildBVH(pivot, end);
    return node;
}