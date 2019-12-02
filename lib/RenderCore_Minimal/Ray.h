#pragma once

class Ray
{
public:
	float3 O;	//ray origin
	float3 E;	//ray direction (normalised)

	Ray::Ray(float3 O, float3 E) : O(O), E(E)
	{
	}
	Ray::~Ray()
	{
	}


};