#pragma once
#include "../gx_vector.h"
#include "gx_material.h"

#include "../binary_tree.h"

class BoundingBox
{
public:
    BoundingBox() = default;
	BoundingBox(const Vector& lbBbox, const Vector& ubBbox);
    
	// Lower bound and upper bound of Bounding Box
	Vector lbBbox, ubBbox;

	bool isInside(const Vector& v) const;
    bool intersect(const Ray& ray) const;
    bool intersect(const Ray& ray, double& t_lower_bound) const;

	void merge(const BoundingBox& other);
	void resizeAddVector(const Vector& v);

    static BoundingBox merge(const BoundingBox& a, const BoundingBox& b);
	static BoundingBox initWithLimits();
};

class ObjectHit;

class Object 
{
public:
	Material material;

    virtual ~Object() {};
    virtual bool intersect(const Ray& ray, ObjectHit& gHit) const = 0;
    virtual bool updateIntersect(const Ray& ray, ObjectHit& gHit) const = 0;

    virtual void transformTranslate(const Vector& delta) = 0;
    virtual void transformScale(double r) = 0;
    virtual void transformRotate(const Vector& axis, double theta) = 0;
    virtual void transformRotate(int axis, double theta) = 0;

protected:
	Object() {}
	Object(Material m) : material(m) {}
};

class ObjectHit 
{
public:
    Vector P;
    Vector N;
    double tHit = std::numeric_limits<double>::max();
	const Object *object_ptr = nullptr;
	
};

class Sphere : public Object
{
public:
    Sphere(const Vector& C, double R);
	Sphere(const Vector& C, double R, const Material& material);

    ~Sphere() = default;

    Vector C;
    double R;
    BoundingBox bbox;

    bool intersect(const Ray& ray, ObjectHit& gHit) const;
	bool updateIntersect(const Ray& ray, ObjectHit& gHit) const;

    void transformTranslate(const Vector& delta);
    void transformScale(double r);
    void transformRotate(const Vector& axis, double theta);
    void transformRotate(int axis, double theta);
private:
	void initBbox();
};

// ----- Triangle meshes ------

class TriangleIndices {
public:
	TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1,
					int nj = -1, int nk = -1, int uvi = -1, int uvj = -1,
					int uvk = -1, int group = -1, bool added = false)
		: vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk),
		ni(ni), nj(nj), nk(nk), group(group){};
		
	int vtxi, vtxj, vtxk; // indices within the vertex coordinates array
	int uvi, uvj, uvk;    // indices within the uv coordinates array
	int ni, nj, nk;       // indices within the normals array
	int group;            // face group
};

struct dataBVH {
	int start; // Index of first triangle
	int end;  // Index of last triangle
	BoundingBox bbox; // Bounding Box
};

class TriangleMesh : public Object{
public:
	TriangleMesh() = delete;
	TriangleMesh(const char *obj);
	TriangleMesh(const char *obj, const Material& material);

	~TriangleMesh() {std::cout << "Triangle Mesh Destructor" << std::endl;}


	void readOBJ(const char *obj);

	static const int BVH_STOP = 5;

	BoundingBox bbox;
	binary_tree::BinaryTree<dataBVH> treeBVH;

	std::vector<TriangleIndices> indices;
	std::vector<Vector> vertices;
	std::vector<Vector> normals;
	std::vector<Vector> uvs;
	std::vector<Vector> vertexcolors;

	bool intersect(const Ray& ray, ObjectHit& gHit) const;
	bool updateIntersect(const Ray& ray, ObjectHit& gHit) const;

	bool fastUpdateIntersect(const Ray& ray, ObjectHit& gHit) const;

	Vector triangleCentroid(int idx) const;

	void transformTranslate(const Vector& delta);
	void transformScale(double r);
	void transformRotate(const Vector& axis, double theta);
	void transformRotate(int axis, double theta);

private:
	// BVH Optimization
	binary_tree::Node<dataBVH>* buildBVH (int start, int end);

	void buildBVH ();
	bool BVHIntersectRoutine(const Ray& ray, ObjectHit& gHit, const binary_tree::Node<dataBVH> *node) const;

	void _verifyTree(binary_tree::Node<dataBVH>* node) const;
	void transformBVH(std::function<void (BoundingBox&)> transform);

};