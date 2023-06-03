#pragma once

#include "../gx_vector.h"

class Material {
public:
    Material() = default;
	Material(const Vector& albedo);

    // White as default color
	Vector albedo = Vector(1., 1. ,1.); 

	bool is_mirror = false;
	bool is_transparent = false;
	double refractive_index = 1.4;

    bool is_light = false;
	double light_intensity = 0;

    void set_light(double light_intensity);
};
