#include "gx_camera.h"

Camera::Camera(const Vector& camera_center, int W, int H, double alpha)
: camera_center(camera_center), W(W), H(H), alpha(alpha)
{
    z_film = -W/(2.*tan(alpha/2.));
}

Ray Camera::ray_to_pixel_center(int i, int j) const
{
    Vector ray_dir = Vector(j - W/2. + 0.5, -i + H/2. + 0.5, z_film);
    return Ray(camera_center, ray_dir.normalized());
}

Ray Camera::ray_to_pixel_gaussian(int i, int j) const
{
    Vector gaussian_vector(0,0,0);
    boxMuller(0.5, gaussian_vector[0], gaussian_vector[1]);

    Vector ray_dir = Vector(j - W/2. + 0.5, -i + H/2. + 0.5, z_film);
    ray_dir += gaussian_vector;

    return Ray(camera_center, ray_dir.normalized());
}

AliasedCamera::AliasedCamera(const Vector& camera_center, int W, int H, double alpha)
: Camera(camera_center, W, H, alpha) 
{}

Ray AliasedCamera::generate_ray(int i, int j) const
{
    return ray_to_pixel_center(i, j);
}

PinholeCamera::PinholeCamera(const Vector& camera_center, int W, int H, double alpha) 
: Camera(camera_center, W, H, alpha) 
{}

Ray PinholeCamera::generate_ray(int i, int j) const
{
    return ray_to_pixel_gaussian(i, j);
}

ProjectiveCamera::ProjectiveCamera(const Vector& camera_center, int W, int H, double alpha,
double focal_distance, double lens_radius)
: Camera(camera_center, W, H, alpha), focal_distance(focal_distance), lens_radius(lens_radius)
{}

Ray ProjectiveCamera::generate_ray(int i, int j) const
{
    Ray center_ray(ray_to_pixel_gaussian(i, j));
    Vector P_focus = camera_center + (focal_distance/std::abs(center_ray.u[2]))*center_ray.u;

    // Generate a random point in the lens - later change this to gx_random
    double r = random_uniform() * sqrt(lens_radius);
    double theta = random_uniform() * 2 * M_PI;

    Vector random_point_lens = camera_center + Vector(r*sin(theta), r*cos(theta), 0.);

    return Ray(random_point_lens, (P_focus-random_point_lens).normalized());
}