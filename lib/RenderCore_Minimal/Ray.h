#pragma once

class Ray
{
public:
	float3 O;	//ray origin
	float3 D;	//ray direction (normalised)

	Ray::Ray(float3 O, float3 D) : O(O), D(D)
	{
	}
	Ray::~Ray()
	{
	}


};