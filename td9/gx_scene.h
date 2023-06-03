#pragma once

#include "../gx_vector.h"
#include "../gx_random.h"
#include "gx_camera.h"
#include "gx_geometry.h"

#include <cassert>

#ifndef EPSILON
    #define EPSILON 1E-3
#endif

static inline double squareDouble(double x) {return x*x; }

static Vector random_cos(const Vector& N){
	// Generate random x, y ,z
	double r1 = random_uniform();
	double r2 = random_uniform();

	double x, y, z;
	x = cos(2*M_PI*r1)*sqrt(1 - r2);
	y = sin(2*M_PI*r1)*sqrt(1 - r2);
	z = sqrt(r2);

	// Generate local frame around N
	Vector T1 = N.orthogonal();
	T1 = T1.normalized();

	Vector T2 = cross(T1, N);
	T2 = T2.normalized();

	return x*T1 + y*T2 + z*N;
}

class Scene {
private:
	bool useIndirectLighting = true;
	std::vector<const Object *> objects;
	
	// Point Light Source
	Vector L;
	double I;

public:
	Scene(const Vector& L, const double& I): L(L), I(I) {};

	void insertObject(const Object* ptr) {
		objects.push_back(ptr);
	}

	void setIndirectLighting(bool newValue) {
		useIndirectLighting = newValue;
	}

	// Physics -- Maybe move to another file later
	double schlick_reflexivity(const Vector& omega_i, const Vector &N, const double &n1, const double& n2){
		// Reflection Coefficient Approximation
		double k0 = ((n1-n2)*(n1-n2))/((n1+n2)*(n1+n2)); 

		return k0 + (1-k0)*std::pow(1 - dot(omega_i, N), 5);
	}

	bool intersect(const Ray& ray, ObjectHit& hitInfo) const
	{	
		hitInfo = ObjectHit();

		for(size_t i = 0; i < objects.size(); i++) {
			ObjectHit possibleHit;
			if (objects[i]->intersect(ray, possibleHit)) {
				if (possibleHit.tHit < hitInfo.tHit) {
					hitInfo = possibleHit;
				}
			}
		}
		return (hitInfo.object_ptr != nullptr);
	}

	bool updateIntersect(const Ray& ray, ObjectHit& hitInfo) const
	{	
		bool hasIntersection = false;
		for (const Object *object_ptr : objects) {
			hasIntersection = object_ptr->updateIntersect(ray, hitInfo) || hasIntersection;
		}
		return hasIntersection;
	}

	bool isVisible(const Vector& P, const Vector& N, const Vector& x) const
	{
		// Returns if point P with normal N can see x directly.
		ObjectHit hitInfo;
		Ray xray = Ray(P+EPSILON*N, (x-P).normalized());

		if (!updateIntersect(xray, hitInfo))
			return true;
		
		return (x-P).norm2() <= squareDouble(hitInfo.tHit);
	}

	Vector directLighting(const ObjectHit& hitInfo) 
	{
		Vector Lo = Vector(0, 0, 0);
		if (!isVisible(hitInfo.P, hitInfo.N, L))
			return Lo;

		Vector omega_i = L - hitInfo.P;
		Vector albedo = hitInfo.object_ptr->material.albedo;

		Lo = (
			(I/(4*M_PI*omega_i.norm2())) *
			(albedo/M_PI) * 
			std::max(0.0, cosine_similarity(hitInfo.N, omega_i))
		);

		return Lo;
	}

	Vector sphereDirectLighting(const ObjectHit& hitInfo, const Sphere *object_ptr)
	{	
		const Material& luminousMaterial = object_ptr->material;

		// Sample point on the light sphere
		const Vector& C = object_ptr->C;
		double R = object_ptr->R;

		Vector D = (hitInfo.P - C).normalized();
		Vector V = random_cos(D);
		Vector x_sphere = C + R * V;

		Vector Lo(0., 0., 0.);

		Vector N_sphere = (x_sphere - C).normalized();
		if (!isVisible(hitInfo.P, hitInfo.N, x_sphere + EPSILON*N_sphere))
			return Lo;

		Vector omega_i = (x_sphere - hitInfo.P).normalized();
		double R2 = (x_sphere - hitInfo.P).norm2();

		double form_factor = (
			std::max(0., dot(omega_i, hitInfo.N)) *
			std::max(0., dot(-omega_i, N_sphere)) /
			(x_sphere - hitInfo.P).norm2()
		);

		double pdf = std::max(0., dot(N_sphere, D))/(M_PI*R2);

		Vector albedo = hitInfo.object_ptr->material.albedo;

		Lo = (
			(albedo/M_PI) *
			(luminousMaterial.light_intensity/(4*M_PI*M_PI*R2)) *
			form_factor /
			pdf
		);

		return Lo;
	}

	Vector indirectLighting(const ObjectHit& hitInfo, int ray_depth){
		if (ray_depth <= 0) return Vector(0,0,0);

		Ray random_ray = Ray(hitInfo.P + EPSILON*hitInfo.N, random_cos(hitInfo.N));
		return hitInfo.object_ptr->material.albedo * getColor(random_ray, ray_depth-1, true);
	}

	Vector getColor(const Ray& ray, int ray_depth, bool last_bounce_diffuse=false){
		if (ray_depth <= 0) return Vector(0,0,0);

		ObjectHit hitInfo;
		if (!updateIntersect(ray, hitInfo))
			return Vector(0.,0.,0.);

		// Hit happens
		Vector& P = hitInfo.P;
		Vector& N = hitInfo.N;
		const Material& materialHit = hitInfo.object_ptr->material;

		if (materialHit.is_light) {
			if (last_bounce_diffuse)
				return Vector(0.,0.,0.);
			else {
				const Sphere *lightSphere_ptr = dynamic_cast<const Sphere *> (hitInfo.object_ptr);
				assert(lightSphere_ptr);

				double R2 = squareDouble(lightSphere_ptr->R);
				Vector light = Vector::initOnesVector() * materialHit.light_intensity;
				return light / (4*M_PI*M_PI*R2);
			}
		}

		if (materialHit.is_mirror) {
			Vector mirror_direction = ray.u - 2 * dot(ray.u, N) * N;
			Ray mirror_ray = Ray(P + EPSILON * N, mirror_direction);
			return getColor(mirror_ray, ray_depth - 1);
		}

		if (materialHit.is_transparent) {
			double n1 = 1; 
			double n2 = materialHit.refractive_index;

			double cos_incidence = dot(ray.u, N);
			bool is_inside_material = (cos_incidence > 0);

			if (is_inside_material) {
				std::swap(n1, n2);
				N *= -1;
				cos_incidence *= -1;
			}
			double n_ratio = n1/n2;
			double rad = 1 - n_ratio*n_ratio*(1 - cos_incidence*cos_incidence);

			if (rad < 0) {
				// total reflection
				Vector mirror_direction = ray.u - 2 * dot(ray.u, N) * N;
				Ray mirror_ray = Ray(P + EPSILON * N, mirror_direction);
				return getColor(mirror_ray, ray_depth - 1);
			}
			Vector transmitted_direction = n_ratio * (ray.u - cos_incidence * N) - sqrt(rad) * N;
			return getColor(Ray(P - EPSILON*N, transmitted_direction), ray_depth-1);
		}

		Vector Lo_direct = directLighting(hitInfo);

		// Check for direct lighting from spheres
		for (auto object_ptr : objects) {
			if (dynamic_cast<const Sphere *> (object_ptr) && object_ptr->material.is_light)
				Lo_direct += sphereDirectLighting(hitInfo, dynamic_cast<const Sphere *> (object_ptr));
		}

		Vector Light = Lo_direct;

		if (useIndirectLighting) {
			Light = Lo_direct + indirectLighting(hitInfo, ray_depth);
		}
		
		return Light;
	}
};