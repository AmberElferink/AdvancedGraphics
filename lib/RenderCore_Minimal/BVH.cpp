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


/*Method that computes the bounding box for a given set of primitives*/
void BVHNode::CalculateBounds( const vector<float4> &vData)
{
	float3 lowest;
	lowest.x = vData[leftFirst].x;
	lowest.y = vData[leftFirst].y;
	lowest.z = vData[leftFirst].z;
	float3 highest;
	highest.x = vData[leftFirst].x;
	highest.y = vData[leftFirst].y;
	highest.z = vData[leftFirst].z;

	int vDataCount = n_objs * 3;
	for ( int i = 1; i < vDataCount; i++ )
	{
		if ( lowest.x > vData[i].x )
			lowest.x = vData[i].x;
		if ( lowest.y > vData[i].y )
			lowest.y = vData[i].y;
		if ( lowest.z > vData[i].z )
			lowest.z = vData[i].z;
		if ( highest.x < vData[i].x )
			highest.x = vData[i].x;
		if ( highest.y < vData[i].y )
			highest.y = vData[i].y;
		if ( highest.z < vData[i].z )
			highest.z = vData[i].z;
	}

	bbox = aabb( lowest, highest );
}

void BVHNode::makeNode(BVH &bvh, unsigned int left_index, unsigned int n_objs)
{


}


//TODO: lege nodes
//TODO: minder dan threshold driehoeken, doe het niet

void BVHNode::SetSplit(BVH &bvh, float &split_point, int &split_index)
{
	float3 D = abs(bbox.bmax3 - bbox.bmin3); //dimension sizes
	if (D.x > D.y)
	{
		if (D.x > D.z) // x is the largest
		{
			quickSortL1byL2X(bvh.indices, bvh.centroids, 0, bvh.nrTriangles - 1);
			split_point = (bbox.bmax3.x - bbox.bmin3.x) / 2;
			split_index = binarySearchX(bvh.indices, bvh.centroids, 0, bvh.nrTriangles - 1, split_point);
		}
		else // z is the largest
		{
			quickSortL1byL2Z(bvh.indices, bvh.centroids, 0, bvh.nrTriangles - 1);
			split_point = (bbox.bmax3.z - bbox.bmin3.z) / 2;
			split_index = binarySearchZ(bvh.indices, bvh.centroids, 0, bvh.nrTriangles - 1, split_point);
		}
	}
	else
	{
		if (D.y > D.z) //y is the largest
		{
			quickSortL1byL2Y(bvh.indices, bvh.centroids, 0, bvh.nrTriangles - 1);
			split_point = (bbox.bmax3.y - bbox.bmin3.y) / 2;
			split_index = binarySearchY(bvh.indices, bvh.centroids, 0, bvh.nrTriangles - 1, split_point);
		}
		else // z is the largest
		{
			quickSortL1byL2Z(bvh.indices, bvh.centroids, 0, bvh.nrTriangles - 1);
			split_point = (bbox.bmax3.z - bbox.bmin3.z) / 2;
			split_index = binarySearchZ(bvh.indices, bvh.centroids, 0, bvh.nrTriangles - 1, split_point);
		}
	}
}

void BVHNode::ConstructRoot(BVH &bvh)
{
	leftFirst = 0;
	n_objs = bvh.nrTriangles;
	(*this).CalculateBounds(bvh.vData);
	
	float splitpoint;
	int split_index;
	SetSplit(bvh, splitpoint, split_index);

	(*this).makeNode(bvh, bvh.poolPtr, split_index);
	bvh.poolPtr++;
	(*this).makeNode(bvh, bvh.poolPtr, bvh.nrTriangles - split_index);
}

void BVH::CalcCentroids()
{
	int nrTriangles = vData.size() / 3;
	centroids.resize(nrTriangles);
	for (int i = 0; i < nrTriangles; i++)
	{
		//v1 v2 and v3 are stored [v1,v2,v3,v1,v2,v3..] so times 3 + 0 give v1, times 3 + 1 gives v2, etc
		float4 centroid = (vData[i * 3] + vData[i * 3 + 1] + vData[i * 3 + 2]) / 3;
		centroids[i] = make_float3(centroid.x, centroid.y, centroid.z);
	}
	quickSortL1byL2Y(indices, centroids, 0, indices.size() - 1);
}

/*Methode that construct a BVH*/
//TODO: build threshold for with no bvh is constructed.
void BVH::ConstructBVH(const float4 *vertexData, const int vertexCount)
{
	vData.resize(vertexCount);
	for (int i = 0; i < vertexCount; i++)
		vData[i] = vertexData[i];
	// create index array
	nrTriangles = vertexCount / 3;
	indices.resize( nrTriangles );
	for ( int i = 0; i < nrTriangles; i++ )
		indices[i] = i;
	// allocate pool for all nodes
	pool.resize( nrTriangles * 2 - 1 );
	root = &pool[0];
	CalcCentroids();
	root->ConstructRoot((*this));
	poolPtr = 2; //start on 2, to keep one node empty next to the root node to align in memory
	// subdivide root node
	left_index = 0;
	right_index = indices.size();
}
