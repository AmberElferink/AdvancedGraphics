#pragma once
#include "common_types.h"
class Camera
{
	Camera(float3 E, float3 V, float d)
	{
		(*this).E = E;
		(*this).V = V;
		(*this).C = E + d * V;
		(*this).P0 = C + (-1, 1, 0);
		(*this).P1 = C + (1, 1, 0);
		(*this).P2 = C + (-1, 1, 0);
	}
	~Camera()
	{

	}
	float3 E = make_float3(0, 0, 0);
	float3 V = make_float3(0, 0, 1);
	float3 C;
	float3 P0;
	float3 P1;
	float3 P2;
};