#include "core_settings.h"

// adapted from Möller–Trumbore intersection algorithm: https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
bool BVHNode::Intersect( const Ray &ray, const CoreTri &triangle, const vector<Material *> &matList, Intersection &intersection )
{
	//TODO: het kan zijn dat een aantal dingen al geprecalculate zijn in CoreTri. Kijk daarnaar voor versnelling
	float3 edge1, edge2, h, s, q;
	float a, f, u, v;
	edge1 = triangle.vertex1 - triangle.vertex0;
	edge2 = triangle.vertex2 - triangle.vertex0;
	h = cross( ray.D, edge2 );
	a = dot( edge1, h );
	if ( a > -EPSILON && a < EPSILON )
		return false; // This ray is parallel to this triangle.
	f = 1.0f / a;
	s = ray.O - triangle.vertex0;
	u = f * dot( s, h );
	if ( u < 0.0f || u > 1.0f )
		return false;
	q = cross( s, edge1 );
	v = f * dot( ray.D, q );
	if ( v < 0.0f || u + v > 1.0f )
		return false;
	// At this stage we can compute t to find out where the intersection point is on the line.
	float t = f * dot( edge2, q );
	if ( t > EPSILON && t < 1.0f / EPSILON ) // ray intersection
	{
		float3 intersectionPoint = ray.O + ray.D * t;
		float3 dir = ray.O - intersectionPoint;
		dir = dir / length( dir );
		float3 normal = make_float3( triangle.Nx, triangle.Ny, triangle.Nz );

		intersection = Intersection( t, intersectionPoint, normal, triangle );
		intersection.material = *matList[triangle.material];
		return true;
	}
	else // This means that there is a line intersection but not a ray intersection.
		return false;
}

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

/* Method that tries to find a split plane for which: A_left * N_left + A_right * N_right < A * N 
   axis = 0 for x-axis, axis = 1 for y-axis and axis = 2 for z-axis*/
void BVHNode::SAH( float &total, int &axis, float &split, const vector<uint> &indices, const vector<aabb> &boundingBoxes, const int leftF )
{
	float current_value = bounds.Area() * count;
	for ( int a = 0; a < 3; a++ ) //do so for every axis
	{
		for ( int i = 0; i < count; i++ )
		{
			split = boundingBoxes[indices[i + leftF]].Center( axis );
			aabb area_left;  //area on left side
			aabb area_right; //area on right side
			uint leftCount = 0;
			uint rightCount = 0;
			bool firstLeft = true;
			bool firstRight = true;

			//check for all elements whether they are to the left or the right of the split plane
			for ( int j = 0; j < count; j++ )
			{
				aabb current_box = boundingBoxes[indices[j + leftF]];
				if ( current_box.Center( axis ) < split ) //case primitive is on the left side
				{
					if ( firstLeft )
					{
						area_left = current_box;
						firstLeft = false;
					}
					else
						area_left.Grow( current_box );
					leftCount++;
				}
				else //case primitive is on the right side
				{
					if ( firstRight )
					{
						area_right = current_box;
						firstRight = false;
					}
					else
						area_right.Grow( current_box );
					rightCount++;
				}
			}
			if ( leftCount < count && rightCount < count )
			{
				//Check whether we make an improvement with the current split or not
				float aleft = area_left.Area();
				float aright = area_right.Area();
				total = (float)leftCount * aleft + (float)rightCount * aright;
				if ( total < current_value )
					return;
			}
		}
		axis = ( axis + 1 ) % 3; //switch axis
	}
	//If no improvement is found, return the current value
}

#define BINNING

#ifdef BINNING
/* Method that divides the three axis in k bins and then compares the 3*(k-1) possible split planes.
   It returns the best split plane based on the SAH-heuristic. 
   based on the article: On fast Construction of SAH-based Bounding Volume Hierarchies, Wald, 2007 */
void BVHNode::Binning( vector<uint> &indices, const vector<aabb> &boundingBoxes, const int leftF, vector<BVHNode> &pool, int &poolPtr )
{
	uint k = 8; //number of bins
	float current_value = bounds.Area() * count; //best SAH-value
	vector<vector<Bin>> bins;
	bins.resize( 3 );
	vector<vector<uint>> bin_indices; //vectors that contains the corresponding bin index for every primitive
	bin_indices.resize( 3 );
	uint best_split = 0; //best split bin
	uint axis; //best split axis

	BVHNode *left;
	left = &pool[poolPtr - 2];
	BVHNode *right;
	right = &pool[poolPtr - 1];

	for ( int a = 0; a < 3; a++ ) //fill bins for every axis
	{
		bins[a].resize( k );
		for ( int i = 0; i < k; i++ )
			bins[a][i].bounds.Reset(); //set aabb of bins to [INF,-INF]
		bin_indices[a].resize( count );
		//precalculated constants
		float k1 = (float)k * ( 1.0f - FLT_EPSILON ) / ( bounds.bmax[a] - bounds.bmin[a] );
		float k0 = bounds.bmin[a];

		//put every primitive in the right bin for the current axis
		for ( int i = 0; i < count; i++ )
		{
			aabb box = boundingBoxes[indices[i + leftF]];
			float p = box.Center( a );
			uint bin_id = k1 * ( p - k0 );
			bins[a][bin_id].Add( box );
			bin_indices[a][i] = bin_id;
		}

		vector<aabb> bbs_right;	//a list with the total bounding boxes of the primitives that are of the right side of the split plane for every possible split
		bbs_right.resize( k - 1 ); //k-1 possible split planes per axis
		bbs_right[k - 2] = ( bins[a][k - 1].bounds );

		//first step is to fill the vector with bounding boxes for the primitives on the right
		for ( int i = k - 3; i > -1; i-- )
		{
			bbs_right[i] = bbs_right[i + 1];
			bbs_right[i].Grow( bins[a][i + 1].bounds );
		}

		aabb bb_left;
		bb_left = bins[a][0].bounds;		//the total bounding box of the primitives that are of the left side of the split plane
		uint count_left = bins[a][0].count; //the total number of primitives that are on the left side of the split plane

		if ( count_left > 0 && count_left < count )
		{
			float sah = bb_left.Area() * count_left + bbs_right[0].Area() * ( count - count_left );
			if ( sah < current_value )
			{
				best_split = 1;
				axis = a;
				left->bounds = bb_left;
				right->bounds = bbs_right[0];
				left->count = count_left;
				right->count = count - count_left;
				current_value = sah;
			}
		}

		for ( int i = 1; i < k - 1; i++ )
		{
			//Add the leftmost bin to the bins on the left
			count_left += bins[a][i].count;
			bb_left.Grow( bins[a][i].bounds );
			if ( count_left > 0 && count_left < count )
			{
				float sah = bb_left.Area() * count_left + bbs_right[i].Area() * ( count - count_left );
				if ( sah < current_value ) //current best solution found
				{
					best_split = i + 1;
					axis = a;
					left->bounds = bb_left;
					right->bounds = bbs_right[i];
					left->count = count_left;
					right->count = count - count_left;
					current_value = sah;
				}
			}
		}
	}

	/* if no improvement is found, we split the current node in half */
	if ( best_split == 0 )
	{
		left->count = ceil( (float)( count / 2 ) );
		left->leftFirst = leftF;
		left->CalculateBounds( boundingBoxes, indices );
		right->count = floor( (float)( count / 2 ) );
		right->leftFirst = leftF + left->count;
		right->CalculateBounds( boundingBoxes, indices );
		left->Subdivide( pool, poolPtr, indices, boundingBoxes );
		right->Subdivide( pool, poolPtr, indices, boundingBoxes );
		return;
	}

	/* swap that is performed with the best_split value */
	uint pointer_low = leftF;
	uint pointer_high = leftF + count - 1;
	uint current = leftF; //the current value in the bin_indices list that we are now considering
	for ( int i = 0; i < count; i++ )
	{
		if (bin_indices[axis][current-leftF] < best_split)
		{
			pointer_low++;
			current = pointer_low;
		}
		else
		{
			uint save = indices[pointer_low];
			indices[pointer_low] = indices[pointer_high];
			indices[pointer_high] = save;
			current = pointer_high;
			pointer_high--;
		}
	}
	left->leftFirst = leftF;
	right->leftFirst = leftF + left->count;
	left->Subdivide( pool, poolPtr, indices, boundingBoxes );
	right->Subdivide( pool, poolPtr, indices, boundingBoxes );
}

/* Method that adds a primitive to a bin and updates the count and aabb accordingly */
void Bin::Add( const aabb &prim )
{
	bounds.Grow( prim );
	count++;
}

#endif

void BVHNode::Partition( vector<uint> &indices, vector<BVHNode> &pool, int &poolPtr, const vector<aabb> &boundingBoxes, const int leftF )
{
	int axis = bounds.LongestAxis();
	float split = 0;
	float total = bounds.Area() * count;

	SAH( total, axis, split, indices, boundingBoxes, leftF );
	if ( total < bounds.Area() * count )
	{
		uint current = leftF;				   //save the index of the first primitive
		uint current_last = leftF + count - 1; //pointer to the last element in the back of the vector that has not been checked
		for ( int i = 0; i < count; i++ )
		{
			if ( boundingBoxes[indices[current]].Center( axis ) < split ) //if left
				current++;
			else //if right, do swap
			{
				//swap elements
				uint save = indices[current];
				indices[current] = indices[current_last];
				indices[current_last] = save;
				current_last--;
			}
		}
		//fill left and right node with values
		BVHNode *left;
		left = &pool[poolPtr - 2];
		left->count = current - leftF;
		left->leftFirst = leftF;
		left->CalculateBounds( boundingBoxes, indices );
		BVHNode *right;
		right = &pool[poolPtr - 1];
		right->count = count - left->count;
		right->leftFirst = leftF + left->count;
		right->CalculateBounds( boundingBoxes, indices );
		left->Subdivide( pool, poolPtr, indices, boundingBoxes );
		right->Subdivide( pool, poolPtr, indices, boundingBoxes );
	}
	else //all the centers lie on the same x,y and z axis, thus we split them randomly
	{
		BVHNode *left;
		left = &pool[poolPtr - 2];
		left->count = ceil( (float)( count / 2 ) );
		left->leftFirst = leftF;
		left->CalculateBounds( boundingBoxes, indices );
		BVHNode *right;
		right = &pool[poolPtr - 1];
		right->count = floor( (float)( count / 2 ) );
		right->leftFirst = leftF + left->count;
		right->CalculateBounds( boundingBoxes, indices );
		left->Subdivide( pool, poolPtr, indices, boundingBoxes );
		right->Subdivide( pool, poolPtr, indices, boundingBoxes );
	}
}

void BVHNode::Subdivide( vector<BVHNode> &pool, int &poolPtr, vector<uint> &indices, const vector<aabb> &boundingBoxes )
{
	if ( count < 3 )
		return; //leaf
	int left = leftFirst;
	leftFirst = poolPtr++; //for the left leaf
	poolPtr++;			   //for the right leaf
	Binning( indices, boundingBoxes, left, pool, poolPtr );
	//Partition( indices, pool, poolPtr, boundingBoxes, left );
}

/*Method that computes the bounding box for a given set of primitives*/
/*TODO:: vertexData omschrijven*/
void BVHNode::CalculateBounds( const vector<aabb> &boundingBoxes, const vector<uint> &indices )
{
	float3 lowest;
	aabb bigBox = boundingBoxes[indices[leftFirst]]; // a big bounding box containing all triangles

	for ( int i = 1; i < count; i++ )
	{
		bigBox.Grow( boundingBoxes[indices[leftFirst + i]] );
	}
	bounds = bigBox;
}

/*Method that traverses trough the nodes of an BVH and returns the closest intersection*/
void BVHNode::Traverse( const Ray &ray, vector<BVHNode> &pool, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection &closest, const vector<Material *> &matList )
{
	if ( !IntersectNode( ray ) )
		return;
	if ( IsLeaf() )
		IntersectPrimitives( ray, indices, triangles, closest, matList );
	else
	{
		pool[leftFirst].Traverse( ray, pool, indices, triangles, closest, matList );
		pool[leftFirst + 1].Traverse( ray, pool, indices, triangles, closest, matList );
	}
}

/*Method that traverses trough the nodes of an BVH and returns the closest intersection*/
bool BVHNode::TraverseToFirst( const Ray &ray, vector<BVHNode> &pool, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection &intersection, const vector<Material *> &matList )
{
	if ( !IntersectNode( ray ) )
		return false;
	if ( IsLeaf() )
	{
		int right = leftFirst + count;
		for ( int i = leftFirst; i < right; i++ )
		{
			if (Intersect(ray, triangles[indices[i]], matList, intersection))
				return true;
		}
		return false;
	}
	else
	{
		return pool[leftFirst].TraverseToFirst( ray, pool, indices, triangles, intersection, matList );
		return pool[leftFirst + 1].TraverseToFirst( ray, pool, indices, triangles, intersection, matList );
	}
}

/* Method that computes the closest intersection for a set of triangles that are contained in one leaf */
void BVHNode::IntersectPrimitives( const Ray &ray, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection &closest, const vector<Material *> &matList )
{
	//compute last index of the array with triangles
	int right = leftFirst + count;
	for ( int i = leftFirst; i < right; i++ )
	{
		Intersection intersection;
		if ( Intersect( ray, triangles[indices[i]], matList, intersection ) )
			if ( intersection.t < closest.t ) //check whether the current intersection is the closest intersection
				closest = intersection;
	}
}

/* Method that checks whether the current ray intersects the bounding box of a given node
   based on: https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection */
bool BVHNode::IntersectNode( const Ray &ray )
{
	float tmin, tmax, tymin, tymax, tzmin, tzmax;

	tmin = ( bounds.MinMax( ray.sign[0], 0 ) - ray.O.x ) * ray.invDir.x;
	tmax = ( bounds.MinMax( 1 - ray.sign[0], 0 ) - ray.O.x ) * ray.invDir.x;
	tymin = ( bounds.MinMax( ray.sign[1], 1 ) - ray.O.y ) * ray.invDir.y;
	tymax = ( bounds.MinMax( 1 - ray.sign[1], 1 ) - ray.O.y ) * ray.invDir.y;

	if ( ( tmin > tymax ) || ( tymin > tmax ) )
		return false;
	if ( tymin > tmin )
		tmin = tymin;
	if ( tymax < tmax )
		tmax = tymax;

	tzmin = ( bounds.MinMax( ray.sign[2], 2 ) - ray.O.z ) * ray.invDir.z;
	tzmax = ( bounds.MinMax( 1 - ray.sign[2], 2 ) - ray.O.z ) * ray.invDir.z;

	if ( ( tmin > tzmax ) || ( tzmin > tmax ) )
		return false;
	if ( tzmin > tmin )
		tmin = tzmin;
	if ( tzmax < tmax )
		tmax = tzmax;

	return true;
}

/*Methode that construct a BVH*/
void BVH::ConstructBVH( const vector<float4> &vertexData, const int vertexCount, const vector<CoreTri> &triangleData )
{
	// create index array
	triangles = triangleData;
	int nrTriangles = vertexCount / 3;
	indices.resize( nrTriangles );
	vector<aabb> boundingBoxes;
	boundingBoxes.resize( nrTriangles );
	for ( int i = 0; i < nrTriangles; i++ )
	{

		aabb boundingBox = aabb( triangleData[i].vertex0, triangleData[i].vertex0 );
		boundingBox.Grow( triangleData[i].vertex1 );
		boundingBox.Grow( triangleData[i].vertex2 );
		float area = boundingBox.Area();
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
	root->CalculateBounds( boundingBoxes, indices );
	root->Subdivide( pool, poolPtr, indices, boundingBoxes );
}
