#pragma once
#include "Ray.h"
#include <iostream>

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
		color = color;
	}
	// data members

	float3 color = make_float3( 0.5 ); // diffuse and reflective material color
	float3 absorption;				   //for dieelectrics
	Texture *texture = 0;			   // texture
	bool metallic = false;
	bool dielectric = false;
	float specularity = 0.8;
	float indexOfRefraction = 1.0003; //air

	float3 GetColor()
	{
		return color;
	}
};

class Intersection
{
  public:
	float t = 10e30;						 // distance from starting point to intersection point
	float3 point = make_float3( 0 );		 // intersection point
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
};

class Bin
{
  public:
	aabb bounds; //bounding box of the primitives in the current bin
	int first;
	int count = 0;
	void Add( const aabb &prim );
};

class BVHNode
{
  public:
	aabb bounds;
	int leftFirst; //left points to the index of the BVHNode in the pool (no leaf) and First is the first index of the primitive contained in the leaf (Leaf)
	int count;	 //number of primitives that are contained in the current node
	inline bool IsLeaf();
	inline int Right();
	void Subdivide( vector<BVHNode> &pool, int &poolPtr, vector<uint> &indices, const vector<aabb> &boundingBoxes );
	void Partition( vector<uint> &indices, vector<BVHNode> &pool, int &poolPtr, const vector<aabb> &boundingBoxes, const int leftF );
	void CalculateBounds( const vector<aabb> &boundingBoxes, const vector<uint> &indices );
	void SAH( float &total, int &axis, float &split, const vector<uint> &indices, const vector<aabb> &boundingBoxes, const int leftF );
	void Binning( vector<uint> &indices, const vector<aabb> &boundingBoxes, const int leftF, vector<BVHNode> &pool, int &poolPtr );
	void Traverse( const Ray &ray, vector<BVHNode> &pool, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection &closest, const vector<Material *> &matList );
	bool TraverseToFirst( const Ray &ray, vector<BVHNode> &pool, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection &closest, const vector<Material *> &matList );
	void IntersectPrimitives( const Ray &ray, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection &closest, const vector<Material *> &matList );
	bool IntersectNode( const Ray &ray );
	static bool Intersect( const Ray &ray, const CoreTri &triangle, const vector<Material *> &matList, Intersection &intersection );
};

class BVH
{
  public:
	//BVH();
	void ConstructBVH( const vector<float4> &vertexData, const int vertexCount, const vector<CoreTri> &triangleData );
	vector<CoreTri> triangles;
	float4 *vertexData;   //list with vertices of all primitives
	vector<uint> indices; //list with indices of all primitives
	vector<BVHNode> pool; //list that contains all the BVH nodes in the BVH
	BVHNode *root;		  //pointer to the root node of the BVH
	int poolPtr;
	//	BVH( const BVH &obj ); // copy constructor
	//	~BVH();
};

//BVH::BVH( const BVH &obj )
//{
//	cout << "Copy constructor allocating root." << endl;
//	root = new BVHNode;
//	*root = *obj.root; // copy the value
//}