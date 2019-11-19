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
	float3 E = make_float3(0, 0, 0); //camera position
	float3 V = make_float3(0, 0, 1); // view direction; must be normalized
	float3 C; // center of screen 
	float3 P0; // upper left corner of screen
	float3 P1; // upper right corner of screen
	float3 P2; // bottom left corner of screen
};