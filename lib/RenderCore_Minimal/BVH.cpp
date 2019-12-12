//#include "core_settings.h"

bool BVHNode::IsLeaf()
{
	if ( count < 3 )
		return true;
	return false;
}

int BVHNode::Right()
{
	return leftFirst + 1;
}

void BVHNode::Partition( vector<uint> &indices, const float4 *primitives, vector<BVHNode> &pool, int &poolPtr, const float4 *vertexData, int leftF )
{
	float split; //x- y- or z-axis split value
	uint current = leftF; //save the index of the first primitive
	uint current_last = leftF + count; 
	for (int i = 0; i < count; i++)
	{
		if ( primitives[indices[current]*3].x < split )
			current++;
		else
		{
			uint save = indices[current];
			indices[current] = indices[current_last];
			indices[current_last] = save;
			current_last--;
		}
	}
	//fill left and right node with values 
	BVHNode *left;
	left = &pool[poolPtr - 1];
	left->count = current - leftFirst;
	left->leftFirst = current;
	left->CalculateBounds( vertexData, left, left->count );
	left->Subdivide( pool, poolPtr, indices, primitives, vertexData );
	BVHNode *right;
	right = &pool[poolPtr];
	right->count = current_last - leftFirst;
	right->leftFirst = count - current; 
	right->CalculateBounds( vertexData, right, right->count );
	right->Subdivide( pool, poolPtr, indices, primitives, vertexData );
}

void BVHNode::Subdivide( vector<BVHNode> &pool, int &poolPtr, vector<uint> &indices, const float4 *primitives, const float4 *vertexData )
{
	if ( count < 3 ) return; //leaf
	int left = leftFirst;
	leftFirst = poolPtr++;
	poolPtr++; //for the right leaf
	Partition(indices,primitives, pool, poolPtr, vertexData, left);
}

/*Method that computes the bounding box for a given set of primitives*/
void BVHNode::CalculateBounds( const float4 *primitives, BVHNode *node, const int count )
{
	float3 lowest;
	lowest.x = primitives[0].x;
	lowest.y = primitives[0].y;
	lowest.z = primitives[0].z;
	float3 highest;
	highest.x = primitives[0].x;
	highest.y = primitives[0].y;
	highest.z = primitives[0].z;

	for ( int i = 1; i < count; i++ )
	{
		if ( lowest.x > primitives[i].x )
			lowest.x = primitives[i].x;
		if ( lowest.y > primitives[i].y )
			lowest.y = primitives[i].y;
		if ( lowest.z > primitives[i].z )
			lowest.z = primitives[i].z;
		if ( highest.x > primitives[i].x )
			highest.x = primitives[i].x;
		if ( highest.y > primitives[i].y )
			highest.y = primitives[i].y;
		if ( highest.z > primitives[i].z )
			highest.z = primitives[i].z;
	}

	node->bounds = aabb( lowest, highest );
}

/*Methode that construct a BVH*/
void BVH::ConstructBVH( const float4 *vertexData, const int vertexCount, const float4 *primitives )
{
	// create index array
	int nrTriangles = vertexCount / 3;
	indices.resize( nrTriangles );
	for ( int i = 0; i < nrTriangles; i++ )
		indices[i] = i;
	// allocate BVH root node
	pool.resize( nrTriangles * 2 - 1 );
	root = &pool[0];
	poolPtr = 2;
	// subdivide root node
	root->leftFirst = 0;
	root->count = nrTriangles;
	root->CalculateBounds( vertexData, root, root->count );
	root->Subdivide( pool, poolPtr, indices, primitives, vertexData );
}
