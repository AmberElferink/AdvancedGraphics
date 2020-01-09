#pragma once

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
/*class Texture
{
  public:
	// constructor / destructor
	Texture() = default;
	Texture( int w, int h ) : width( w ), height( h ) { pixels = (uint *)MALLOC64( w * h * sizeof( uint ) ); }
	~Texture() { FREE64( pixels ); }
	// data members
	int width = 0, height = 0;
	uint *pixels = 0;
};*/

/*// -----------------------------------------------------------
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
		color = color;
	}
	// data members

	float3 color = make_float3( 0.5 ); // diffuse and reflective material color
	float3 absorption; //for dieelectrics
	Texture *texture = 0;			   // texture
	bool metallic = false;
	bool dielectric = false;
	float specularity = 0.8;
	float indexOfRefraction = 1.0003; //air

	float3 GetColor()
	{
		return color;
	}
};*/

//  +-----------------------------------------------------------------------------+
//  |  Mesh                                                                       |
//  |  Minimalistic mesh storage.                                           LH2'19|
//  +-----------------------------------------------------------------------------+
class Mesh
{
  public:

	vector<float4> vertices;   // vertex data received via SetGeometry
	int vcount = 0;			   // vertex count
	vector<CoreTri> triangles; // 'fat' triangle data
	BVH bvh;
	//Mesh( const Mesh &obj ); // copy constructor
};

class Light
{
  public:
	Light() = default;
	// data members
	bool pointLight = false;
	bool directionalLight = false;
	bool spotLight = false;
	bool areaLight = false;
	float3 direction = make_float3( 0, 0, 0 );
	float3 position = make_float3( 0, 0, 0 );
	float3 radiance = make_float3( 0, 0, 0 );
	float cosInner = 0;
	float cosOuter = 0;
	CoreLightTri triangle;
	float energy = 1.0;
	int dummy = 1;
};

/*class Intersection
{
  public:
	float t;						 // distance from starting point to intersection point
	float3 point;							 // intersection point
	float3 norm = make_float3( -1, -1, -1 ); // normal at intersection point
	Material material;
	CoreTri triangle;

	Intersection( const float t, const float3 &point, const float3 &norm, const CoreTri &triangle ) : t( t ), point( point ), norm( norm )
	{
		//material = triangle.material;
	}
	Intersection()
	{
	}
}*/

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
	vector<BVH> bvh;
	BitmapFloat *buffer;
	Scene scene;
	int2 probePos; //mouse position for which you want to get more info

	// constructor / destructor
	Raytracer() = default;
	// methods
	void Init();
	void Reinit( int w, int h, Surface *screen );
	void Render( const mat4 &transform );
	//bool Intersect( const Ray &ray, const CoreTri &triangle, Intersection &intersection );
	bool IsOccluded( const Ray &ray, const Light &light );
	bool viewLight( const Intersection &intersection, const Light &light, float3 &lightVector );
	bool viewDirLight( const Intersection &intersection, const Light &light, float3 &lightVector );
	int viewSpotLight( const Intersection &intersection, const Light &light, float3 &lightVector );
	bool viewAreaLight( const Intersection &intersection, Light &light );
	float3 randomPointTri( const CoreLightTri &triangle );
	void rayTrace( Bitmap *screen, const ViewPyramid &view, const int targetTextureID );
	void rayTraceBlock( const ViewPyramid &view, Bitmap *screen, const int targetTextureID, int lineStart, int lineEnd );
	void rayTraceLine( Bitmap *screen, const ViewPyramid &view, const int targetTextureID, const int lineNr );
	void rayTraceBlockAVX(const ViewPyramid &view, Bitmap *screen, const int targetTextureID, int lineStart, int lineEnd);
	void rayTraceLineAVX(Bitmap *screen, const ViewPyramid &view, const int targetTextureID, const int lineNr);
	void rayTraceInPackets(Bitmap *screen, const ViewPyramid &view, const int targetTextureID, const int lineNr);
	uint FloatToIntColor( float3 floatColor );
	Intersection nearestIntersection( const Ray &ray );
	void nearestIntersection(const Ray8 &ray, Intersection8 &closest);
	void nearestIntersection( Rays &r, Intersections &closests);
	Bitmap *rayTraceRandom( const ViewPyramid &view, const int targetTextureID, int &frameCounter );
	void TextureColor( Intersection &intersection, const CoreTri &triangle, uint &color );


	//Pathtracer methodes
	float3 Sample(const Ray &ray);
	void pathTrace(Bitmap *screen, const ViewPyramid &view, const int targetTextureID, uint sampleCount);

	float3 Reflect( const Ray &ray, const Intersection &intersection, int reflectionDepth );
	float3 calcDielectric( Ray ray, Intersection intersection, const Intersection prevIntersection, int reflectionDepth, const float n1 = 1.0002f ); //only adjust n1 if previous trace is also a dielectric material
	float Fresnel( const float cosi, const float ncalc, const float n1, const float n2 );
	float3 DirectIllumination( Intersection intersection );
	float3 DiffuseReflection( float3 N );
	float3 Trace( const Ray &ray, const Intersection prevIntersection, int reflectionDepth ); //default: air
	Color8 Trace(Ray8 &ray, const Intersection8 prevIntersection, int reflectionDepth); //default: air
	void Trace(Rays &ray, Indices I, const Intersections prevIntersection, int reflectionDepth); //default: air
	void storeBVH();
};
} // namespace lh2core