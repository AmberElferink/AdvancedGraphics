#pragma once
#include "Ray.h"
//#include "common_types.h"

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
	~Surface() { FREE64( pixels ); /* assuming we used MALLOC64 to create the buffer */ }
	int width = 0, height = 0;
	uint *pixels = 0;
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
	Texture( int w, int h ) : width( w ), height( h ) { pixels = (uint *)MALLOC64( w * h * sizeof( uint ) ); }
	~Texture() { FREE64( pixels ); }
	// data members
	int width = 0, height = 0;
	uint *pixels = 0;
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
	Material( float3 color )
	{
		diffuse = color;
	}
	// data members

	float3 diffuse = make_float3( 1 ); // diffuse material color
	Texture *texture = 0;			   // texture

	float3 GetColor()
	{
		return diffuse;
	}
};

//  +-----------------------------------------------------------------------------+
//  |  Mesh                                                                       |
//  |  Minimalistic mesh storage.                                           LH2'19|
//  +-----------------------------------------------------------------------------+
class Mesh
{
  public:
	float4 *vertices = 0;   // vertex data received via SetGeometry
	int vcount = 0;			// vertex count
	CoreTri *triangles = 0; // 'fat' triangle data
};

class Light
{
  public:
	Light() = default;
	// data members
	float3 position = make_float3( 0, 0, 0 );
	float3 radiance = make_float3( 0, 0, 0 );
	float energy = 1.0;
	int dummy = 1;
};

class Intersection
{
  public:
	float t;								 // distance from starting point to intersection point
	float3 point;							 // intersection point
	float3 norm = make_float3( -1, -1, -1 ); // normal at intersection point
	Material *material;

	Intersection( const float t, const float3 &point, const float3 &norm, const CoreTri &triangle ) : t( t ), point( point ), norm( norm )
	{
		//material = triangle.material;
	}
	Intersection()
	{
	}
};

class Scene
{
  public:
	vector<Material *> matList;
	vector<Texture *> texList;
	vector<Light> lightList;
	vector<Mesh> meshList;
};

// -----------------------------------------------------------
// Rasterizer class
// rasterizer
// implements a basic, but fast & accurate software rasterizer
// -----------------------------------------------------------
class Raytracer
{
  public:
	Scene scene;
	// constructor / destructor
	Raytracer() = default;
	// methods
	void Init();
	void Reinit( int w, int h, Surface *screen );
	void Render( const mat4 &transform );
	bool Intersect( const Ray &ray, const CoreTri &triangle, Intersection &intersection );
	bool IsOccluded( const Ray &ray, const Light &light );
	bool viewLight( float3 I, const Light &light, float3 lightVector );
	void rayTrace( Bitmap *screen, const ViewPyramid &view );
	uint FloatToIntColor( float3 floatColor );
};
} // namespace lh2core