#include "core_settings.h"

//bool BVHNode::isLeaf()
//{
//	if ( n_objs < 3 )
//		return true;
//	return false;
//}
//
//int BVHNode::right()
//{
//	return index + 1;
//}
//
//void BVHNode::Partition(vector<uint> &indices, const uint start, const uint end)
//{
//	
//}
//
//void BVHNode::Subdivide( vector<BVHNode> &pool, int &poolPtr, vector<uint> &indices )
//{
//	if ( count < 3 ) return; //leaf
//	leftFirst = poolPtr++;
//	poolPtr++; //update poolPtr for the right node
//	Partition(indices,start,end);
//	Subdivide( pool, poolPtr, indices ); //left
//	poolPtr--; //reset for the right child
//	Subdivide( pool, poolPtr, indices ); //right
//}

///*Method that computes the bounding box for a given set of primitives*/
//void BVHNode::CalculateBounds( const float4 *primitives, BVHNode *node, const int count )
//{
//	float3 lowest;
//	lowest.x = primitives[0].x;
//	lowest.y = primitives[0].y;
//	lowest.z = primitives[0].z;
//	float3 highest;
//	highest.x = primitives[0].x;
//	highest.y = primitives[0].y;
//	highest.z = primitives[0].z;
//
//	for ( int i = 1; i < count; i++ )
//	{
//		if ( lowest.x > primitives[i].x )
//			lowest.x = primitives[i].x;
//		if ( lowest.y > primitives[i].y )
//			lowest.y = primitives[i].y;
//		if ( lowest.z > primitives[i].z )
//			lowest.z = primitives[i].z;
//		if ( highest.x > primitives[i].x )
//			highest.x = primitives[i].x;
//		if ( highest.y > primitives[i].y )
//			highest.y = primitives[i].y;
//		if ( highest.z > primitives[i].z )
//			highest.z = primitives[i].z;
//	}
//
//	node->bounds = aabb( lowest, highest );
//}

/*Methode that construct a BVH*/
void BVH::ConstructBVH(const float4 *vertexData, const int vertexCount)
{
	triangles t;
	t.vData = vertexData;
	int w = 0;
	printf("hoi");
	//// create index array
	//int nrTriangles = vertexCount / 3;
	//indices.resize( nrTriangles );
	//for ( int i = 0; i < nrTriangles; i++ )
	//	indices[i] = i;
	//// allocate BVH root node
	//pool.resize( nrTriangles * 2 - 1 );
	//root = &pool[0];
	//poolPtr = 2;
	//// subdivide root node
	//root->index = 0;
	//root->n_objs = nrTriangles;
	//root->CalculateBounds( vertexData, root, root->n_objs );
	//root->Subdivide( pool, poolPtr, indices );
}
