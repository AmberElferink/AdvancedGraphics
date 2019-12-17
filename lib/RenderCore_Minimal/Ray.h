#pragma once

class Ray
{
public:
	float3 O;	//ray origin
	float3 D;	//ray direction (normalised)
	float3 I = make_float3(1); //Intensity of the light remaining after traveling for every color, this diminishes by absorption in a dielectric.
	float3 invDir;
	int sign[3]; 

	Ray::Ray(float3 O, float3 D) : O(O), D(D)
	{
		invDir = 1 / D;
		sign[0] = ( invDir.x < 0 );
		sign[1] = ( invDir.y < 0 );
		sign[2] = ( invDir.z < 0 );
	}
	Ray::~Ray()
	{
	}


};