#pragma once

class Ray
{
public:
	float3 O;	//ray origin
	float3 E;	//ray direction
	float d = -1;	//distance to the nearest primitive

	Ray::Ray(float3 O, float3 E) : O(O), E(E)
	{
	}
	Ray::~Ray()
	{
	}


};