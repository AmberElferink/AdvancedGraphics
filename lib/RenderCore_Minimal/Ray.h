#pragma once
#define COHERENTTRAVERSAL //packet traversal, 4 rays at once
#define RAYPACKETSIZE 4 //8 simd rays per packet, so 32 rays total
#define TOTALPACKETSIZE 8 * RAYPACKETSIZE 
static union { __m256 trueMask8; float trueMask[8]; };
class Ray
{
public:
	float3 O;	//ray origin
	float3 D;	//ray direction (normalised)
	float3 I = make_float3(1); //Intensity of the light remaining after traveling for every color, this diminishes by absorption in a dielectric.
	float3 color = make_float3(0.0f); //rayColor up to the current bounce
	float3 T = make_float3(1.0f);
	float3 recDir; //reciprocal is used multiple times per ray
	int signX; //this is also used multiple times per ray
	int signY;
	int signZ;

	Ray::Ray(float3 O, float3 D) : O(O), D(D)
	{
		recDir = 1 / D;
		signX = ( recDir.x < 0 );
		signY = ( recDir.y < 0 );
		signZ = ( recDir.z < 0 );
	}
	inline Ray Ray::Reflect(const float3 &N, const float3 &point) const
	{
		float3 reflectedDir = D - 2 * dot(D, N) * N;
		return  Ray(point + 2 * EPSILON * reflectedDir, reflectedDir);
	}
	Ray::~Ray()
	{
	}
};
#ifdef COHERENTTRAVERSAL

class Color8
{
public:
	union { __m256 r8; float b[8]; };
	union { __m256 g8; float g[8]; };
	union { __m256 b8; float r[8]; };
	Color8(float r, float g, float b)
	{
		r8 = _mm256_set1_ps(r);
		g8 = _mm256_set1_ps(g);
		b8 = _mm256_set1_ps(b);
	}
	Color8()
	{
		r8 = _mm256_setzero_ps();
		g8 = _mm256_setzero_ps();
		b8 = _mm256_setzero_ps();
	}
};

class Colors
{
public:
	Color8 colors[RAYPACKETSIZE];
};

class ALIGN(32) Ray8 //Eight rays in one, AVX for SIMD stuff. 
{
public:

	union { __m256 ox8; float ox[8]; };
	union { __m256 oy8; float oy[8]; };
	union { __m256 oz8; float oz[8]; };

	union { __m256 dx8; float dx[8]; };
	union { __m256 dy8; float dy[8]; };
	union { __m256 dz8; float dz[8]; };

	union { __m256 recDirX8; float recDirX[8]; };
	union { __m256 recDirY8; float recDirY[8]; };
	union { __m256 recDirZ8; float recDirZ[8]; };

	union { __m256 signX8; float signX[8]; };
	union { __m256 signY8; float signY[8]; };
	union { __m256 signZ8; float signZ[8]; };

	union { __m256 activeMask8; float deadMask[8]; };

	Color8 color;

	Ray8::Ray8(const float3* O, const float3* D) //8 origins and 8 directions
	{
		color = Color8(1.f, 0.f, 1.f);
		for (int i = 0; i < 8; i++)
		{
			ox[i] = O[i].x;
			oy[i] = O[i].y;
			oz[i] = O[i].z;

			dx[i] = D[i].x;
			dy[i] = D[i].y;
			dz[i] = D[i].z;
		}

			recDirX8 = _mm256_rcp_ps(dx8); // 1/dx
			recDirY8 = _mm256_rcp_ps(dy8);
			recDirZ8 = _mm256_rcp_ps(dz8);

			signX8 = _mm256_cmp_ps(recDirX8, _mm256_setzero_ps(), _CMP_LT_OS);
			signY8 = _mm256_cmp_ps(recDirY8, _mm256_setzero_ps(), _CMP_LT_OS);
			signZ8 = _mm256_cmp_ps(recDirZ8, _mm256_setzero_ps(), _CMP_LT_OS);

			//somehow, just throwing the truemask set in raytracer constructor or rendercore init doesn't work correctly.
			trueMask8 = _mm256_cmp_ps(_mm256_setzero_ps(), _mm256_setzero_ps(), _CMP_EQ_OS);
			activeMask8 = trueMask8;// _mm256_cmp_ps(_mm256_setzero_ps(), _mm256_setzero_ps(), _CMP_EQ_OS);// trueMask8;
	}
	Ray8::Ray8()
	{
		trueMask8 = _mm256_cmp_ps(_mm256_setzero_ps(), _mm256_setzero_ps(), _CMP_EQ_OS);
		activeMask8 = trueMask8;
	}
	//only put true, the bool itself will not be checked
	Ray8::Ray8(bool raysStartDead)
	{
		activeMask8 = _mm256_setzero_ps();
	}
	Ray8::~Ray8()
	{
	}
};

//group of Ray8 simd rays for packet traversal
class Rays
{
public:
	Ray8 rays[RAYPACKETSIZE];
	
	Rays() = default;
	//only put true, the bool itself will not be checked
	Rays(bool raysStartDead)
	{
		for (int i = 0; i < RAYPACKETSIZE; i++)
		{
			rays[i].activeMask8 = _mm256_setzero_ps();
		}
	}
	//int ia = RAYPACKETSIZE; //one past last active ray (rays at and behind rays[I[ia]] do not intersect)
};

//this class is a wrapper around the indices corresponding to rays. This array is sorted so the active rays are in the front and the inactive rays are in the back.
class Indices
{
public:
	int I[RAYPACKETSIZE];
	Indices()
	{
		for (int i = 0; i < RAYPACKETSIZE; i++)
		{
			I[i] = i;
		}
	}
};


#endif