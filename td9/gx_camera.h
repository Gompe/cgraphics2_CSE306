#pragma once

#include "../gx_vector.h"
#include "../gx_random.h"

class Camera
{
public:
    Camera(const Vector& camera_center, int W, int H, double alpha);
    Vector camera_center;
    int W;
    int H;
    double alpha;
    // z-distance between camera center and film
    double z_film; 

    virtual Ray generate_ray(int i, int j) const = 0;
    
    Ray ray_to_pixel_center(int i, int j) const;
    Ray ray_to_pixel_gaussian(int i, int j) const;
};

class AliasedCamera : public Camera
{
public:
    AliasedCamera(const Vector& camera_center, int W, int H, double alpha);
    Ray generate_ray(int i, int j) const;
};

class PinholeCamera : public Camera
{
public:
    PinholeCamera(const Vector& camera_center, int W, int H, double alpha);
    Ray generate_ray(int i, int j) const;
};

class ProjectiveCamera : public Camera
{
public:
    ProjectiveCamera(const Vector& camera_center, int W, int H, double alpha, 
            double focal_distance, double radius_aperture);
    
    double focal_distance;
    double lens_radius;

    Ray generate_ray(int i, int j) const;
};