#pragma once

class Ray
{
public:
	float3 O;	//ray origin
	float3 D;	//ray direction (normalised)
	float3 I = make_float3(1); //Intensity of the light remaining after traveling for every color, this diminishes by absorption in a dielectric.

	Ray::Ray(float3 O, float3 D) : O(O), D(D)
	{
	}
	Ray::~Ray()
	{
	}


};