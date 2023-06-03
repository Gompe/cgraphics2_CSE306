#include "gx_material.h"

Material::Material(const Vector& albedo)
: albedo(albedo)
{}

void Material::set_light(double light_intensity)
{
    is_light = (light_intensity != 0.);
	this->light_intensity = light_intensity;
}