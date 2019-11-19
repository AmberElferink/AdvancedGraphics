#pragma once
#include "common_types.h"

class Ray
{
	Ray(float3 O, float3 E, float d) : O(O), E(E), d(d)
	{
	}
	~Ray()
	{
	}
	float3 O;
	float3 E;
	float d;
};