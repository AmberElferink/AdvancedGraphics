#pragma once
#include "Ray.h"
//#include "common_types.h"


class Intersection
{
  public:
	float t;								 // distance from starting point to intersection point
	float3 point;							 // intersection point
	float3 norm = make_float3( -1, -1, -1 ); // normal at intersection point
	int material = 0x0000ff;				 //red for now, figure this out later. Probably a property somewhere in the triangle but can't find where

	Intersection( const float t, const float3 &point, const float3 &norm, const CoreTri &triangle ) : t( t ), point( point ), norm( norm )
	{
		//material = triangle.material;
	}
	Intersection()
	{
	}
};

// adapted from Möller–Trumbore intersection algorithm: https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
bool Intersect( const Ray &ray, const CoreTri &triangle, Intersection &intersection )
{
	//TODO: het kan zijn dat een aantal dingen al geprecalculate zijn in CoreTri. Kijk daarnaar voor versnelling
	float3 vertex0 = triangle.vertex0;
	float3 vertex1 = triangle.vertex1;
	float3 vertex2 = triangle.vertex2;
	float3 edge1, edge2, h, s, q;
	float a, f, u, v;
	edge1 = vertex1 - vertex0;
	edge2 = vertex2 - vertex0;
	h = cross( ray.E, edge2 );
	a = dot( edge1, h );
	if ( a > -EPSILON && a < EPSILON )
		return false; // This ray is parallel to this triangle.
	f = 1.0 / a;
	s = ray.O - vertex0;
	u = f * dot( s, h );
	if ( u < 0.0 || u > 1.0 )
		return false;
	q = cross( s, edge1 );
	v = f * dot( ray.E, q );
	if ( v < 0.0 || u + v > 1.0 )
		return false;
	// At this stage we can compute t to find out where the intersection point is on the line.
	float t = f * dot( edge2, q );
	if ( t > EPSILON && t < 1 / EPSILON ) // ray intersection
	{
		float3 intersectionPoint = ray.O + ray.E * t;
		float3 normal = make_float3( triangle.Nx, triangle.Ny, triangle.Nz );
		intersection = Intersection( t, intersectionPoint, normal, triangle );
		return true;
	}
	else // This means that there is a line intersection but not a ray intersection.
		return false;
}

// The rest of this file is adapted from rasterizer

namespace lh2core
{

	// -----------------------------------------------------------
	// Surface class
	// bare minimum
	// -----------------------------------------------------------
	class Surface
	{
	public:
		Surface() = default;
		~Surface() { FREE64(pixels); /* assuming we used MALLOC64 to create the buffer */ }
		int width = 0, height = 0;
		uint* pixels = 0;
	};

	// -----------------------------------------------------------
	// Texture class
	// encapsulates a palettized pixel surface with pre-scaled
	// palettes for fast shading
	// -----------------------------------------------------------
	class Texture
	{
	public:
		// constructor / destructor
		Texture() = default;
		Texture(int w, int h) : width(w), height(h) { pixels = (uint*)MALLOC64(w * h * sizeof(uint)); }
		~Texture() { FREE64(pixels); }
		// data members
		int width = 0, height = 0;
		uint* pixels = 0;
	};

	// -----------------------------------------------------------
	// Material class
	// basic material properties
	// -----------------------------------------------------------
	class Material
	{
	public:
		// constructor / destructor
		Material() = default;
		// data members
		uint diffuse = 0xffffffff;		// diffuse material color
		Texture* texture = 0;			// texture
	};

	// -----------------------------------------------------------
	// SGNode class
	// scene graph node, with convenience functions for translate
	// and transform; base class for Mesh
	// -----------------------------------------------------------
	class SGNode
	{
	public:
		enum { SG_TRANSFORM = 0, SG_MESH };
		// constructor / destructor
		~SGNode()
		{
			for (int s = (int)child.size(), i = 0; i < s; i++)
			{
				for (int j = i + 1; j < s; j++) if (child[j] == child[i]) child[j] = 0;
				delete child[i];
			}
		}
		// methods
		void SetPosition(float3& pos) { mat4& M = localTransform; M[3] = pos.x, M[7] = pos.y, M[11] = pos.z; }
		float3 GetPosition() { mat4& M = localTransform; return make_float3(M[3], M[7], M[11]); }
		void Render(const mat4& transform);
		virtual int GetType() { return SG_TRANSFORM; }
		// data members
		mat4 localTransform;
		vector<SGNode*> child;
	};


	// -----------------------------------------------------------
	// Scene class
	// owner of the scene graph;
	// owner of the material and texture list
	// -----------------------------------------------------------
	class Scene
	{
	public:
		// constructor / destructor
		Scene() = default;
		~Scene();
		// data members

		SGNode* root = 0;
		vector<Material*> matList;
		vector<Texture*> texList;
	};

	// -----------------------------------------------------------
	// Rasterizer class
	// rasterizer
	// implements a basic, but fast & accurate software rasterizer
	// -----------------------------------------------------------
	class Raytracer
	{
		public:
			// constructor / destructor
			Raytracer() = default;
			// methods
			void Init();
			void Reinit(int w, int h, Surface* screen);
			void Render(const mat4& transform);
			// data members
			Scene scene;
			static float* zbuffer;
			static float4 frustum[5];
	};
}