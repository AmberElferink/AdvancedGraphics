#pragma once
#include "Ray.h"
#include <iostream>

#ifdef COHERENTTRAVERSAL
const __m256 EPS8 = _mm256_set1_ps(EPSILON);
const __m256 MINUSEPS8 = _mm256_set1_ps(-EPSILON);
const __m256 ONE8 = _mm256_set1_ps(1.0f);
#endif

extern int rayNr;
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
	// data members

	float3 color = make_float3( 0.5 ); // diffuse and reflective material color
	float3 absorption;				   //for dieelectrics
	Texture *texture = 0;			   // texture
	bool metallic = false;
	bool dielectric = false;
	float specularity = 0.8;
	float indexOfRefraction = 1.0003; //air

		// constructor / destructor /copy constructor
	Material() = default;
	Material(float3 color)
	{
		color = color;
	}
	Material(const Material &mat)
	{
		color = mat.color;
		absorption = mat.absorption;
		texture = mat.texture;
		metallic = mat.metallic;
		dielectric = mat.dielectric;
		specularity = mat.specularity;
		indexOfRefraction = mat.indexOfRefraction;
	}

	float3 GetColor()
	{
		return color;
	}
};

//frustrum for culling triangles out of the four corner rays in packet traversal.
class Frustrum
{
private:
	// the normal and offsets for all the side planes of the frustrum
	float nx[4];
	float ny[4];
	float nz[4];
	float d[4]; //offset plane
public:
	Frustrum(const Rays &r)
	{
		// 0,0     plane0        0,7  
		// ____________________
		// |ray0               | ray1
		// | plane3            |
		// |                   | plane1
		// |  Packetsize-1,0   |
		//______________________ PACKETSIZE -1, 7
		//        plane2
		//calculate the normals of the planes around the frustrum
		// ni = di cross d(i+1) %4
		cross(r.rays[0].dx[0], r.rays[0].dy[0], r.rays[0].dz[0], 
			r.rays[0].dx[7], r.rays[0].dy[7], r.rays[0].dz[7], 
			nx[0], ny[0], nz[0]); //output
		cross( r.rays[0].dx[7], r.rays[0].dy[7], r.rays[0].dz[7],
			r.rays[RAYPACKETSIZE - 1].dx[7], r.rays[RAYPACKETSIZE - 1].dy[7], r.rays[RAYPACKETSIZE - 1].dz[7],
		    nx[1], ny[1], nz[1]);
		cross(r.rays[RAYPACKETSIZE - 1].dx[7], r.rays[RAYPACKETSIZE - 1].dy[7], r.rays[RAYPACKETSIZE - 1].dz[7],
			r.rays[RAYPACKETSIZE - 1].dx[0], r.rays[RAYPACKETSIZE - 1].dy[0], r.rays[RAYPACKETSIZE - 1].dz[0],
		    nx[2], ny[2], nz[2]);
		cross(r.rays[RAYPACKETSIZE - 1].dx[0], r.rays[RAYPACKETSIZE - 1].dy[0], r.rays[RAYPACKETSIZE - 1].dz[0],
			r.rays[0].dx[0], r.rays[0].dy[0], r.rays[0].dz[0],
			nx[3], ny[3], nz[3]);

		//calculate the offsets of the planes around the frustrum
		// di = oi dot ni
		// dot definitie: float c = return a.x * b.x + a.y * b.y + a.z * b.z;
		d[0] = r.rays[0].ox[0] * nx[0] + r.rays[0].oy[0] * ny[0] + r.rays[0].oz[0] * nz[0];
		d[1] = r.rays[0].ox[7] * nx[1] + r.rays[0].oy[7] * ny[1] + r.rays[0].oz[7] * nz[1];
		d[2] = r.rays[RAYPACKETSIZE -1].ox[7] * nx[2] + r.rays[RAYPACKETSIZE - 1].oy[7] * ny[2] + r.rays[RAYPACKETSIZE - 1].oz[7] * nz[2];
		d[3] = r.rays[RAYPACKETSIZE - 1].ox[0] * nx[3] + r.rays[RAYPACKETSIZE - 1].oy[0] * ny[3] + r.rays[RAYPACKETSIZE - 1].oz[0] * nz[3];
	}

	bool Includes(const aabb &bbox) const
	{
		//all boundingbox corners must be checked:
		// 0: bmin
 		// 1: bmin xy bmax z
		// 2. bmin y  bmax xz
		// 3. bmin yz bmax x
		// 4. bmin xz bmax y
		// 5. bmin x  bmax yz
		// 6. bmax
		// 7. bmin z  bmax xy
		// Paper A. Reshetov, A. Soupikov, and J. Hurley. Multi-level ray tracing algorithm. ACM TOG SIGGRAPH 05, 24(3), 2005.
		// this paper contains the knowledge to do this in SSE (no pseudo), but for now the easy way
		float H[8]; //if H > 0 for a corner, the corner is outside the frustrum
		for(int i = 0; i < 4; i++) //loop over all planes in the frustrum
		{
			// Hi = ni dot pk - bi

			//Check for bbox corners 0 to 7 as stated above:
			H[0] = (nx[i] * bbox.bmin3.x + ny[i] * bbox.bmin3.y + nz[i] * bbox.bmin3.z) - d[i];
			H[1] = (nx[i] * bbox.bmin3.x + ny[i] * bbox.bmin3.y + nz[i] * bbox.bmax3.z) - d[i];
			H[2] = (nx[i] * bbox.bmax3.x + ny[i] * bbox.bmin3.y + nz[i] * bbox.bmax3.z) - d[i];
			H[3] = (nx[i] * bbox.bmax3.x + ny[i] * bbox.bmin3.y + nz[i] * bbox.bmin3.z) - d[i];
			H[4] = (nx[i] * bbox.bmin3.x + ny[i] * bbox.bmax3.y + nz[i] * bbox.bmin3.z) - d[i];
			H[5] = (nx[i] * bbox.bmin3.x + ny[i] * bbox.bmax3.y + nz[i] * bbox.bmax3.z) - d[i];
			H[6] = (nx[i] * bbox.bmax3.x + ny[i] * bbox.bmax3.y + nz[i] * bbox.bmax3.z) - d[i];
			H[7] = (nx[i] * bbox.bmax3.x + ny[i] * bbox.bmax3.y + nz[i] * bbox.bmin3.z) - d[i];
			//if Hi > 0, it is outside the plane. If all pk (corners) are outside the same plane, this plane culls the polyhedron.
			for (int i = 0; i < 8; i++)
			{
				if (H[i] <= 0)
					return true; //if only one point of the bbox is in the frustrum, the frustrum includes it.
			}
			//if none of the points gave true, it is culled by this plane the frustrum.
			return false;
		}

	}
};

class Intersection
{
  public:
	float t = 10e30;								 // distance from starting point to intersection point
	float3 point = make_float3(0);							 // intersection point
	float3 norm = make_float3( -1, -1, -1 ); // normal at intersection point
	Material material;
	CoreTri triangle;

    Intersection( const float t, const float3 &point, const float3 &norm, const CoreTri &triangle ) : t( t ), point( point ), norm( norm )
	{
		//material = triangle.material;
	}
	Intersection::Intersection(const Intersection& inter)
	{
		t = inter.t;
		point = inter.point;
		norm = inter.norm;
		material = inter.material;
		triangle = inter.triangle;
	}
	Intersection()
	{
	}
};

class Intersection8
{
public:
	Intersection intersections[8];
	union { float t[8]; __m256 t8; };
	Intersection8(float* t8)
	{
		for (int i = 0; i < 8; i++)
		{
			t[i] = t8[i];
		}
	}
	Intersection8()
	{
		for (int i = 0; i < 8; i++)
		{
			t[i] = 10e30;
		}
	}

	void SetMaterialsZero()
	{
		for (int i = 0; i < 8; i++)
		{
			intersections[i].material = make_float3(0);
		}
	}
};
//cluster of intersections for packet traversal
class Intersections
{
public:
	Intersection8 inter[RAYPACKETSIZE];
};
#ifdef COHERENTTRAVERSAL
class ALIGN(32) BVHNode
{
public:
	aabb bounds;
	int leftFirst; //left points to the index of the BVHNode in the pool (no leaf) and First is the first index of the primitive contained in the leaf (Leaf)
	uint count;	 //number of primitives that are contained in the current node
	inline bool IsLeaf();
	inline int Right();
	void Subdivide(vector<BVHNode> &pool, int &poolPtr, vector<uint> &indices, const vector<aabb> &boundingBoxes);
	void Partition(vector<uint> &indices, vector<BVHNode> &pool, int &poolPtr, const vector<aabb> &boundingBoxes, const int leftF);
	void CalculateBounds(const vector<aabb> &boundingBoxes, const vector<uint> &indices);
	void SAH(float &total, int &axis, float &split, const vector<uint> &indices, const vector<aabb> &boundingBoxes, const int leftF);
	void Traverse(const Ray &ray, vector<BVHNode> &pool, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection &closest, const vector<Material *> &matList);
	void Traverse(const Ray8 &rays, vector<BVHNode> &pool, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection8 &closest, const vector<Material *> &matList);
	void Traverse(Rays &r, int ia, Indices I, vector<BVHNode> &pool, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersections &closests, const vector<Material *> &matList);
	void Traverse(Rays &r, const Frustrum &fr, int ia, Indices I, vector<BVHNode> &pool, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersections &closests, const vector<Material *> &matList);
	bool TraverseToFirst(const Ray &ray, vector<BVHNode> &pool, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection &intersection, const vector<Material *> &matList);
	void IntersectPrimitives(const Ray &ray, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection &closest, const vector<Material *> &matList);
	void IntersectPrimitives(const Ray8 &rays, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection8 &closest, const vector<Material *> &matList);
	void IntersectPrimitives(const Rays &r, int ia, const Indices &I, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersections &closests, const vector<Material *> &matList);
	bool IntersectNode(const Ray &ray);
	bool IntersectNode(const Ray8 &rays);
	int partRays(Rays &r, int ia, Indices &I);
	int partRays(Rays &r, const Frustrum &fr, int ia, Indices &I);
	static bool Intersect(const Ray &ray, const CoreTri &triangle, const vector<Material*> &matList, Intersection &intersection);
	static bool Intersect(const Ray8 &ray, const CoreTri &triangle, const vector<Material*> &matList, Intersection8 &intr);
	static bool IntersectClosest(const Ray8 &ray, const CoreTri &triangle, const vector<Material*> &matList, Intersection8 &intr);

};
#else
class ALIGN(32) BVHNode
{
  public:
	aabb bounds;
	int leftFirst; //left points to the index of the BVHNode in the pool (no leaf) and First is the first index of the primitive contained in the leaf (Leaf)
	uint count;	 //number of primitives that are contained in the current node
	inline bool IsLeaf();
	inline int Right();
	void Subdivide( vector<BVHNode> &pool, int &poolPtr, vector<uint> &indices, const vector<aabb> &boundingBoxes );
	void Partition( vector<uint> &indices, vector<BVHNode> &pool, int &poolPtr, const vector<aabb> &boundingBoxes, const int leftF );
	void CalculateBounds( const vector<aabb> &boundingBoxes, const vector<uint> &indices );
	void SAH( float &total, int &axis, float &split, const vector<uint> &indices, const vector<aabb> &boundingBoxes, const int leftF );
	void Traverse( const Ray &ray, vector<BVHNode> &pool, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection &closest, const vector<Material *> &matList );
	void Traverse( const Ray8 &rays, vector<BVHNode> &pool, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection8 &closest, const vector<Material *> &matList );
	bool TraverseToFirst(const Ray &ray, vector<BVHNode> &pool, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection &closest, const vector<Material *> &matList);
	void IntersectPrimitives( const Ray &ray, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection &closest, const vector<Material *> &matList );
	void IntersectPrimitives(const Ray8 &rays, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection8 &closest, const vector<Material *> &matList);
	bool IntersectNode( const Ray &ray );
	bool IntersectNode( const Ray8 &rays );
	static bool Intersect( const Ray &ray, const CoreTri &triangle, const vector<Material*> &matList, Intersection &intersection );
	static bool Intersect(const Ray8 &ray, const CoreTri &triangle, const vector<Material*> &matList, Intersection8 &intr);
};
#endif



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