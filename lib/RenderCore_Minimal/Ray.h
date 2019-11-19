#pragma once
#include "common_types.h"

class Ray
{
public:
	Ray(float3 O, float3 E, float d) : O(O), E(E), d(d)
	{
	}
	~Ray()
	{
	}

	float3 O;	//ray origin
	float3 E;	//ray direction
	float d;	//distance to the nearest primitive
};