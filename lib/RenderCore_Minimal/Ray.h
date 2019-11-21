#pragma once
#include "common_types.h"

class Ray
{
public:
	Ray(float3 O, float3 E) : O(O), E(E)
	{
	}
	~Ray()
	{
	}

	float3 O;	//ray origin
	float3 E;	//ray direction
	float d = -1;	//distance to the nearest primitive
};