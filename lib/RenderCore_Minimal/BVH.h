#pragma once

class BVHNode
{
  public:
	aabb bounds;
	int leftFirst; //left points to the index of the BHNode in the pool (no leaf) and First is the first index of the primitive contained in the leaf (Leaf)
	int count; //number of primitives that are contained in the current node
	inline bool IsLeaf();
	inline int Right();
	void Subdivide( vector<BVHNode> &pool, int &poolPtr, vector<uint> &indices );
	void Partition();
	void CalculateBounds( const float4 *primitives, BVHNode *node, const int count );
};

class BVH
{
	void ConstructBVH( const float4 *triangles, const int vertexCount );
	float4 *vertexData;		  //list with vertices of all primitives
	vector<float3> centroids; //list with centroids of all primitives
	vector<uint> indices;	 //list with indices of all primitives
	vector<BVHNode> pool;     //list that contains all the BVH nodes in the BVH
	BVHNode *root;             //pointer to the root node of the BVH
	int poolPtr;
};