#pragma once
class Triangle
{
public:
	float4 v1;
	float4 v2;
	float4 v3;
};

class BVH; //forward declaration: there will be a class BVH

class BVHNode
{
private:
	aabb bbox;
	uint leftFirst;  // if leaf == true, index to left child node.
				// else if leaf == true, index to first triangle in the triangle vector.
	uint n_objs; //number of primitives that are contained in the current node
	inline uint right();
	void CalculateBounds(const vector<float4> &vData);
	void SetSplit(BVH &bvh, float &split_point, int &split_index);
public:
	inline bool isLeaf();
	void makeLeaf(unsigned int index_, unsigned int n_objs);
	void makeNode(BVH &bvh, unsigned int left_index, unsigned int n_objs);
	//n_bjs in makenode is voor debug only
	void setAABB(aabb &bbox_);

	uint getIndex() { return leftFirst; }
	uint getNObjs() { return n_objs; }
	aabb &getAABB() { return bbox; }
	void ConstructRoot(BVH &bvh);
};

class BVH
{
public:
	void ConstructBVH(const float4 *triangles, const int vertexCount);
private:
	friend class BVHNode; //BVHNode can access private BVH parameters

	void CalcCentroids();
	void ConstructRecursive(int left_index, int right_index, aabb box, BVHNode *node, int depth);
	vector<float3> centroids; //list with centroids of all primitives
	vector<uint> indices;	 //list with indices of all primitives
	vector<BVHNode> pool;     //list that contains all the BVH nodes in the BVH
	vector<float4> vData; 
	BVHNode *root;             //pointer to the root node of the BVH
	int poolPtr;
	int nrTriangles;
	int left_index;
	int right_index;
	
};
