#include "core_settings.h"

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

void BVHNode::Partition( vector<uint> &indices, vector<BVHNode> &pool, int &poolPtr, const vector<aabb> boundingBoxes, int leftF )
{
	int longest_side = 0;
	int length_x = abs( bounds.bmax3.x - bounds.bmin3.x );
	int length_y = abs( bounds.bmax3.y - bounds.bmin3.y );
	int length_z = abs( bounds.bmax3.z - bounds.bmin3.z );
	if ( length_y > length_x && length_y > length_z )
		longest_side = 1;
	else if ( length_z > length_x )
		longest_side = 2;

	uint current = leftF; //save the index of the first primitive
	uint current_last = leftF + count - 1;

	if (longest_side == 0)
	{
		float split = ( bounds.bmin3.x + bounds.bmax3.x ) / 2; //x- y- or z-axis split value
		for ( int i = 0; i < count; i++ )
		{
			if ( boundingBoxes[indices[current]].Center( 0 ) < split )
				current++;
			else
			{
				//swap elements
				uint save = indices[current];
				indices[current] = indices[current_last];
				indices[current_last] = save;
				current_last--;
			}
		}
	}
	else if (longest_side == 1)
	{
		float split = ( bounds.bmin3.y + bounds.bmax3.y ) / 2; //x- y- or z-axis split value
		for ( int i = 0; i < count; i++ )
		{
			if ( boundingBoxes[indices[current]].Center( 1 ) < split )
			current++;
			else
			{
			//swap elements
			uint save = indices[current];
			indices[current] = indices[current_last];
			indices[current_last] = save;
			current_last--;
			}
		}
	}
	else
	{
		float split = ( bounds.bmin3.z + bounds.bmax3.z ) / 2; //x- y- or z-axis split value
		for ( int i = 0; i < count; i++ )
		{
			if ( boundingBoxes[indices[current]].Center( 2 ) < split )
				current++;
			else
			{
				//swap elements
				uint save = indices[current];
				indices[current] = indices[current_last];
				indices[current_last] = save;
				current_last--;
			}
		}
	}


	//fill left and right node with values
	BVHNode *left;
	left = &pool[poolPtr - 2];
	left->count = current - leftF;
	left->leftFirst = leftF;
	left->CalculateBounds( boundingBoxes, left->count, indices );
	BVHNode *right;
	right = &pool[poolPtr - 1];
	right->count = count - left->count;
	right->leftFirst = leftF + left->count;
	right->CalculateBounds( boundingBoxes, right->count, indices );
	left->Subdivide( pool, poolPtr, indices, boundingBoxes );
	right->Subdivide( pool, poolPtr, indices, boundingBoxes );
}

void BVHNode::Subdivide( vector<BVHNode> &pool, int &poolPtr, vector<uint> &indices, const vector<aabb> boundingBoxes )
{
	if ( count <= 4 ) 
		return; //leaf
	int left = leftFirst;
	leftFirst = poolPtr++;
	poolPtr++; //for the right leaf
	Partition( indices, pool, poolPtr, boundingBoxes, left );
}

/*Method that computes the bounding box for a given set of primitives*/
/*TODO:: vertexData omschrijven*/
void BVHNode::CalculateBounds( const vector<aabb> boundingBoxes, const int count, vector<uint> indices )
{
	float3 lowest;
	aabb bigBox = boundingBoxes[0]; // a big bounding box containing all triangles

	for ( int i = 1; i < count; i++ )
	{
		bigBox.Grow( boundingBoxes[indices[leftFirst + i]] );
	}
	bounds = bigBox;
}

/*Methode that construct a BVH*/
void BVH::ConstructBVH( const vector<float4> vertexData, const int vertexCount )
{
	// create index array
	int nrTriangles = vertexCount / 3;
	indices.resize( nrTriangles );
	vector<aabb> boundingBoxes;
	boundingBoxes.resize( nrTriangles );
	for ( int i = 0; i < nrTriangles; i++ )
	{
		aabb boundingBox = aabb( vertexData[3*i], vertexData[3*i+1] );
		boundingBox.Grow( vertexData[3 * i + 2] );
		indices[i] = i;
		boundingBoxes[i] = boundingBox;
	}
	// allocate BVH root node
	pool.resize( nrTriangles * 2 - 1 );
	root = &pool[0];
	poolPtr = 2;
	// subdivide root node
	root->leftFirst = 0;
	root->count = nrTriangles;
	root->CalculateBounds( boundingBoxes, root->count, indices );
	root->Subdivide( pool, poolPtr, indices, boundingBoxes );
}
