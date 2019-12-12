#pragma once
class Triangle
{
	float4 v1;
	float4 v2;
	float4 v3;
};

class BVHNode
{
private:
	aabb bbox;
	uint index;  // if leaf == true, index to left child node.
				// else if leaf == true, index to first triangle in the triangle vector.
	uint n_objs; //number of primitives that are contained in the current node
	inline uint right();
	void CalculateBounds(const float4 *primitives, BVHNode *node, const int count);
public:
	inline bool isLeaf();
	void makeLeaf(unsigned int index_, unsigned int n_objs);
	void makeNode(unsigned int left_index, unsigned int n_objs);
	//n_bjs in makenode is voor debug only
	void setAABB(aabb &bbox_);

	uint getIndex() { return index; }
	uint getNObjs() { return n_objs; }
	aabb &getAABB() { return bbox; }

};

class BVH
{
public:

	void ConstructBVH(const float4 *triangles, const int vertexCount);
	void ConstructRecursive(int left_index, int right_index, aabb box, BVHNode *node, int depth);
	vector<float3> centroids; //list with centroids of all primitives
	vector<uint> indices;	 //list with indices of all primitives
	vector<BVHNode> pool;     //list that contains all the BVH nodes in the BVH
	BVHNode *root;             //pointer to the root node of the BVH
	int poolPtr;
	union triangles { const float4 *vData; const Triangle *triangleList; };
};