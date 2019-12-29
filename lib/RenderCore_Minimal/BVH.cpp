#include "core_settings.h"



const __m256 true256 = _mm256_set1_ps(0xFFFFFFFF);
// adapted from Möller–Trumbore intersection algorithm: https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
bool BVHNode::Intersect( const Ray &ray, const CoreTri &triangle, const vector<Material*> &matList, Intersection &intersection )
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
		float3 normal = make_float3( triangle.Nx, triangle.Ny, triangle.Nz );

		intersection = Intersection( t, intersectionPoint, normal, triangle );
		intersection.material = *matList[triangle.material];
		return true;
	}
	else // This means that there is a line intersection but not a ray intersection.
		return false;
}

// adapted from Möller–Trumbore intersection algorithm: https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
// adapted SSE code to AVX from slides Jacco Bikker "SIMD recap" 2016/2017
bool BVHNode::Intersect(const Ray8 &ray8, const CoreTri &triangle, const vector<Material*> &matList, Intersection8 &intr)
{
	__m256 e1x8 = _mm256_set1_ps(triangle.vertex1.x - triangle.vertex0.x);
	__m256 e1y8 = _mm256_set1_ps(triangle.vertex1.y - triangle.vertex0.y);
	__m256 e1z8 = _mm256_set1_ps(triangle.vertex1.z - triangle.vertex0.z);
	__m256 e2x8 = _mm256_set1_ps(triangle.vertex2.x - triangle.vertex0.x);
	__m256 e2y8 = _mm256_set1_ps(triangle.vertex2.y - triangle.vertex0.y);
	__m256 e2z8 = _mm256_set1_ps(triangle.vertex2.z - triangle.vertex0.z);
	__m256 Px8 = _mm256_sub_ps(_mm256_mul_ps(ray8.dy8, e2z8), _mm256_mul_ps(ray8.dz8, e2y8));
	__m256 Py8 = _mm256_sub_ps(_mm256_mul_ps(ray8.dz8, e2x8), _mm256_mul_ps(ray8.dx8, e2z8));
	__m256 Pz8 = _mm256_sub_ps(_mm256_mul_ps(ray8.dx8, e2y8), _mm256_mul_ps(ray8.dy8, e2x8));
	__m256 det8 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(e1x8, Px8), _mm256_mul_ps(e1y8, Py8)), _mm256_mul_ps(e1z8, Pz8));
	// det <= -EPS|| det >= EPS
	__m256 mask1 = _mm256_or_ps(_mm256_cmp_ps(det8, MINUSEPS8, _CMP_LE_OQ), _mm256_cmp_ps(det8, EPS8, _CMP_GE_OQ));

	__m256 inv_det8 = _mm256_rcp_ps(det8);
	__m256 Tx8 = _mm256_sub_ps(ray8.ox8, _mm256_set1_ps(triangle.vertex0.x));
	__m256 Ty8 = _mm256_sub_ps(ray8.oy8, _mm256_set1_ps(triangle.vertex0.y));
	__m256 Tz8 = _mm256_sub_ps(ray8.oz8, _mm256_set1_ps(triangle.vertex0.z));
	__m256 u8 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(Tx8, Px8), _mm256_mul_ps(Ty8, Py8)), _mm256_mul_ps(Tz8, Pz8)), inv_det8);
	__m256 mask2 = _mm256_and_ps(_mm256_cmp_ps(u8, _mm256_setzero_ps(), _CMP_GE_OQ), _mm256_cmp_ps(u8, ONE8, _CMP_LE_OQ));
	__m256 Qx8 = _mm256_sub_ps(_mm256_mul_ps(Ty8, e1z8), _mm256_mul_ps(Tz8, e1y8));
	__m256 Qy8 = _mm256_sub_ps(_mm256_mul_ps(Tz8, e1x8), _mm256_mul_ps(Tx8, e1z8));
	__m256 Qz8 = _mm256_sub_ps(_mm256_mul_ps(Tx8, e1y8), _mm256_mul_ps(Ty8, e1x8));
	__m256 v8 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(ray8.dx8, Qx8), _mm256_mul_ps(ray8.dy8, Qy8)), _mm256_mul_ps(ray8.dz8, Qz8)), inv_det8);
	__m256 mask3 = _mm256_and_ps(_mm256_cmp_ps(v8, _mm256_setzero_ps(), _CMP_GE_OQ), _mm256_cmp_ps(_mm256_add_ps(u8, v8), ONE8, _CMP_LE_OQ));

	union { __m256 t8; float t[8]; };
	t8 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(e2x8, Qx8), _mm256_mul_ps(e2y8, Qy8)), _mm256_mul_ps(e2z8, Qz8)), inv_det8);
	__m256 mask8 = _mm256_cmp_ps(t8, _mm256_setzero_ps(), _CMP_GE_OQ);
	//__m256 mask5 = _mm256_cmp_ps(t8, rays.t8, _CMP_LT_OQ); test if distance is shorter than to previous noted intersection. Happens in nearestIntersection
	__m256 combined = _mm256_and_ps(_mm256_and_ps(_mm256_and_ps(mask1, mask2), mask3), mask8);
	//ray8.t8 = _mm_blendv_ps(ray8.t8, t8, combined); store result in ray for each alive ray.
	if (_mm256_movemask_ps(combined) == 0)
		return false;
	else
	{
		//we currently want to give back intersections for all alive rays. in a masked array
		for (int i = 0; i < 8; i++)
		{
			float3 intersectionPoint = make_float3(ray8.ox[i] + ray8.dx[i] * t[i], ray8.oy[i] + ray8.dy[i] * t[i], ray8.oz[i] + ray8.dz[i] * t[i]);
			float3 normal = make_float3(triangle.Nx, triangle.Ny, triangle.Nz);

			intr.intersections[i] = Intersection(t[i], intersectionPoint, normal, triangle);
			intr.intersections[i].material = *matList[triangle.material];
			int w = 0;
		}
		return true;
	}

}

bool BVHNode::IntersectClosest(const Ray8 &ray8, const CoreTri &triangle, const vector<Material*> &matList, Intersection8 &closest)
{
	__m256 e1x8 = _mm256_set1_ps(triangle.vertex1.x - triangle.vertex0.x);
	__m256 e1y8 = _mm256_set1_ps(triangle.vertex1.y - triangle.vertex0.y);
	__m256 e1z8 = _mm256_set1_ps(triangle.vertex1.z - triangle.vertex0.z);
	__m256 e2x8 = _mm256_set1_ps(triangle.vertex2.x - triangle.vertex0.x);
	__m256 e2y8 = _mm256_set1_ps(triangle.vertex2.y - triangle.vertex0.y);
	__m256 e2z8 = _mm256_set1_ps(triangle.vertex2.z - triangle.vertex0.z);
	__m256 Px8 = _mm256_sub_ps(_mm256_mul_ps(ray8.dy8, e2z8), _mm256_mul_ps(ray8.dz8, e2y8));
	__m256 Py8 = _mm256_sub_ps(_mm256_mul_ps(ray8.dz8, e2x8), _mm256_mul_ps(ray8.dx8, e2z8));
	__m256 Pz8 = _mm256_sub_ps(_mm256_mul_ps(ray8.dx8, e2y8), _mm256_mul_ps(ray8.dy8, e2x8));
	__m256 det8 = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(e1x8, Px8), _mm256_mul_ps(e1y8, Py8)), _mm256_mul_ps(e1z8, Pz8));
	// det <= -EPS|| det >= EPS
	__m256 mask1 = _mm256_or_ps(_mm256_cmp_ps(det8, MINUSEPS8, _CMP_LE_OQ), _mm256_cmp_ps(det8, EPS8, _CMP_GE_OQ));

	__m256 inv_det8 = _mm256_rcp_ps(det8);
	__m256 Tx8 = _mm256_sub_ps(ray8.ox8, _mm256_set1_ps(triangle.vertex0.x));
	__m256 Ty8 = _mm256_sub_ps(ray8.oy8, _mm256_set1_ps(triangle.vertex0.y));
	__m256 Tz8 = _mm256_sub_ps(ray8.oz8, _mm256_set1_ps(triangle.vertex0.z));
	__m256 u8 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(Tx8, Px8), _mm256_mul_ps(Ty8, Py8)), _mm256_mul_ps(Tz8, Pz8)), inv_det8);
	__m256 mask2 = _mm256_and_ps(_mm256_cmp_ps(u8, _mm256_setzero_ps(), _CMP_GE_OQ), _mm256_cmp_ps(u8, ONE8, _CMP_LE_OQ));
	__m256 Qx8 = _mm256_sub_ps(_mm256_mul_ps(Ty8, e1z8), _mm256_mul_ps(Tz8, e1y8));
	__m256 Qy8 = _mm256_sub_ps(_mm256_mul_ps(Tz8, e1x8), _mm256_mul_ps(Tx8, e1z8));
	__m256 Qz8 = _mm256_sub_ps(_mm256_mul_ps(Tx8, e1y8), _mm256_mul_ps(Ty8, e1x8));
	__m256 v8 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(ray8.dx8, Qx8), _mm256_mul_ps(ray8.dy8, Qy8)), _mm256_mul_ps(ray8.dz8, Qz8)), inv_det8);
	__m256 mask3 = _mm256_and_ps(_mm256_cmp_ps(v8, _mm256_setzero_ps(), _CMP_GE_OQ), _mm256_cmp_ps(_mm256_add_ps(u8, v8), ONE8, _CMP_LE_OQ));
	union { __m256 t8; float t[8];};
	t8 = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(e2x8, Qx8), _mm256_mul_ps(e2y8, Qy8)), _mm256_mul_ps(e2z8, Qz8)), inv_det8);
	__m256 mask8 = _mm256_cmp_ps(t8, _mm256_setzero_ps(), _CMP_GE_OQ);
	//__m256 mask5 = _mm256_cmp_ps(t8, closest.t8, _CMP_LT_OQ); //test if distance is shorter than to previous noted intersection. Happens in nearestIntersection
	union {__m256 combined8; float combined[8];};
	combined8 = _mm256_and_ps(_mm256_and_ps(_mm256_and_ps(mask1, mask2), mask3), mask8);//_mm256_and_ps(_mm256_and_ps(_mm256_and_ps(_mm256_and_ps(mask1, mask2), mask3), mask8), mask5);
	//ray8.t8 = _mm_blendv_ps(ray8.t8, t8, combined); store result in ray for each alive ray.
	
	int finalMask = _mm256_movemask_ps(combined8);
	if (finalMask == 0)
		return false;
	else
	{
		closest.t8 = _mm256_blendv_ps(closest.t8, t8, combined8);
		//we currently want to give back intersections for all alive rays. in a masked array
		for (int i = 0; i < 8; i++)
		{
			if (isnan(abs(combined[i])) )
			{
				float3 intersectionPoint = make_float3(ray8.ox[i] + ray8.dx[i] * t[i], ray8.oy[i] + ray8.dy[i] * t[i], ray8.oz[i] + ray8.dz[i] * t[i]);
				float3 normal = make_float3(triangle.Nx, triangle.Ny, triangle.Nz);
				closest.intersections[i] = Intersection(t[i], intersectionPoint, normal, triangle);
				closest.intersections[i].material = *matList[triangle.material];
			}

		}
		return true;
	}

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
		for ( uint i = 0; i < count; i++ )
		{
			split = boundingBoxes[indices[i + leftF]].Center(axis);
			aabb area_left;  //area on left side
			aabb area_right; //area on right side
			uint leftCount = 0;
			uint rightCount = 0;
			bool firstLeft = true;
			bool firstRight = true;

			//check for all elements whether they are to the left or the right of the split plane
			for ( uint j = 0; j < count; j++ )
			{
				aabb current_box = boundingBoxes[indices[j + leftF]];
				if ( current_box.Center( axis ) < split ) //case primitive is on the left side
				{
					if (firstLeft)
					{
						area_left = current_box;
						firstLeft = false;
					}
					else area_left.Grow( current_box );
					leftCount++;
				}
				else //case primitive is on the right side
				{
					if (firstRight)
					{
						area_right = current_box;
						firstRight = false;
					}
					else area_right.Grow( current_box );
					rightCount++;
				}
			}
			if ( leftCount < count && rightCount < count )
			{
			//Check whether we make an improvement with the current split or not
			float aleft = area_left.Area();
		    float aright = area_right.Area();
			total = (float)leftCount * aleft + (float)rightCount * aright;
			if ( total < current_value)
				return;
			}
		}
		axis = ( axis + 1 ) % 3; //switch axis
	}
	//If no improvement is found, return the current value
}

void BVHNode::Partition( vector<uint> &indices, vector<BVHNode> &pool, int &poolPtr, const vector<aabb> &boundingBoxes, const int leftF )
{
	int axis = bounds.LongestAxis();
	float split = 0;
	float total = bounds.Area() * count;

	SAH( total, axis, split, indices, boundingBoxes, leftF );
	if ( total < bounds.Area() * count )
	{
		uint current = leftF; //save the index of the first primitive
		uint current_last = leftF + count - 1; //pointer to the last element in the back of the vector that has not been checked
		for ( uint i = 0; i < count; i++ )
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
		left->count = ceil((float)(count / 2));
		left->leftFirst = leftF;
		left->CalculateBounds( boundingBoxes, indices );
		BVHNode *right;
		right = &pool[poolPtr - 1];
		right->count = floor( (float)(count / 2));
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
	poolPtr++; //for the right leaf
	Partition( indices, pool, poolPtr, boundingBoxes, left );
}

/*Method that computes the bounding box for a given set of primitives*/
/*TODO:: vertexData omschrijven*/
void BVHNode::CalculateBounds( const vector<aabb> &boundingBoxes, const vector<uint> &indices )
{
	aabb bigBox = boundingBoxes[indices[leftFirst]]; // a big bounding box containing all triangles

	for ( uint i = 1; i < count; i++ )
	{
		bigBox.Grow( boundingBoxes[indices[leftFirst + i]] );
	}
	bounds = bigBox;
}

/*Method that traverses trough the nodes of an BVH and returns the closest intersection*/
void BVHNode::Traverse(const Ray &ray, vector<BVHNode> &pool, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection &closest, const vector<Material *> &matList)
{
	if (!IntersectNode(ray))
		return;
	if (IsLeaf())
		IntersectPrimitives(ray, indices, triangles, closest, matList);
	else
	{
		pool[leftFirst].Traverse(ray, pool, indices, triangles, closest, matList);
		pool[leftFirst + 1].Traverse(ray, pool, indices, triangles, closest, matList);
	}
}
/*Method that traverses trough the nodes of an BVH and returns the closest intersection*/
void BVHNode::Traverse( const Ray8 &ray8, vector<BVHNode> &pool, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection8 &closest, const vector<Material *> &matList )
{
	if ( !IntersectNode( ray8 ) )
		return;
	if ( IsLeaf() )
		IntersectPrimitives( ray8, indices, triangles, closest, matList );
	else
	{
		pool[leftFirst].Traverse( ray8, pool, indices, triangles, closest, matList);
		pool[leftFirst + 1].Traverse( ray8, pool, indices, triangles, closest, matList );
	}
}

/*Method that traverses trough the nodes of an BVH and returns if it intersected and the first intersection, useful for shadowrays*/
bool BVHNode::TraverseToFirst(const Ray &ray, vector<BVHNode> &pool, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection &intersection, const vector<Material *> &matList)
{
	if (!IntersectNode(ray))
		return false;
	if (IsLeaf())
	{
		int right = leftFirst + count;
		for (int i = leftFirst; i < right; i++)
		{
			if (Intersect(ray, triangles[indices[i]], matList, intersection))
				return true;
		}
		return false;
	}
	else
	{
		pool[leftFirst].Traverse(ray, pool, indices, triangles, intersection, matList);
		pool[leftFirst + 1].Traverse(ray, pool, indices, triangles, intersection, matList);
	}
}

/* Method that computes the closest intersection for a set of triangles that are contained in one leaf */
void BVHNode::IntersectPrimitives( const Ray &ray, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection &closest, const vector<Material *> &matList )
{
    //compute last index of the array with triangles
	int right = leftFirst + count;
	for (int i = leftFirst; i < right; i++)
	{
		Intersection intersection;
		if(Intersect( ray, triangles[indices[i]], matList, intersection ))
			if ( intersection.t < closest.t ) //check whether the current intersection is the closest intersection
				closest = intersection;
	}	
}

/* Method that computes the closest intersection for a set of triangles that are contained in one leaf */
void BVHNode::IntersectPrimitives(const Ray8 &rays, const vector<uint> &indices, const vector<CoreTri> &triangles, Intersection8 &closest, const vector<Material *> &matList)
{
	//compute last index of the array with triangles
	int right = leftFirst + count;
	for (int i = leftFirst; i < right; i++)
	{
		IntersectClosest(rays, triangles[indices[i]], matList, closest);
	}
}

/* Method that checks whether the current ray intersects the bounding box of a given node
   based on: https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection */
bool BVHNode::IntersectNode(const Ray &ray)
{
	float tmin, tmax, tymin, tymax, tzmin, tzmax;

	tmin = ( bounds.MinMax(ray.signX,0) - ray.O.x ) * ray.recDir.x;
	tmax = ( bounds.MinMax(1 - ray.signX,0) - ray.O.x ) * ray.recDir.x;
	tymin = ( bounds.MinMax(ray.signY,1) - ray.O.y ) * ray.recDir.y;
	tymax = ( bounds.MinMax(1 - ray.signY,1) - ray.O.y ) * ray.recDir.y;

	if ( ( tmin > tymax ) || ( tymin > tmax ) )
		return false;
	if ( tymin > tmin )
		tmin = tymin;
	if ( tymax < tmax )
		tmax = tymax;

	tzmin = ( bounds.MinMax(ray.signZ, 2) - ray.O.z ) * ray.recDir.z;
	tzmax = ( bounds.MinMax(1 - ray.signZ, 2) - ray.O.z ) * ray.recDir.z;

	if ( ( tmin > tzmax ) || ( tzmin > tmax ) )
		return false;

	return true;
}



/* Method that checks whether the current ray intersects the bounding box of a given node
   based on: https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection 
   translated to AVX*/

bool BVHNode::IntersectNode(const Ray8 &ray8)
{
	__m256 tmin, tmax, tymin, tymax, tzmin, tzmax;

	tmin = _mm256_mul_ps(_mm256_sub_ps(bounds.MinMax( /*flip signX8*/  ray8.signX8, 0),           ray8.ox8), ray8.recDirX8);
	tmax = _mm256_sub_ps(_mm256_sub_ps(bounds.MinMax( _mm256_andnot_ps(ray8.signX8, true256), 0), ray8.ox8), ray8.recDirX8);
	tymin = _mm256_sub_ps(_mm256_sub_ps(bounds.MinMax(                 ray8.signY8, 1),			 ray8.oy8), ray8.recDirY8);
	tymax = _mm256_sub_ps(_mm256_sub_ps(bounds.MinMax(_mm256_andnot_ps(ray8.signY8, true256), 1), ray8.oy8), ray8.recDirY8);

	__m256 tminGTtymax = _mm256_cmp_ps(tmin, tymax, _CMP_GT_OS);
	__m256 tyminGTtmax = _mm256_cmp_ps(tymin, tmax, _CMP_GT_OS);
	//or the above two. If any of them is > 0 (use movemask according to slides), continue, otherwise return false
	int checkY = _mm256_movemask_ps(_mm256_or_ps(tminGTtymax, tyminGTtmax));
	if (checkY == 0)
		return false;

	tmin = _mm256_max_ps(tymin, tmin);
	tmax = _mm256_min_ps(tymax, tmax);
	
	tzmin = _mm256_sub_ps(_mm256_sub_ps(bounds.MinMax(                 ray8.signZ8, 2),			 ray8.oz8), ray8.recDirZ8);
	tzmax = _mm256_sub_ps(_mm256_sub_ps(bounds.MinMax(_mm256_andnot_ps(ray8.signZ8, true256), 2), ray8.oz8), ray8.recDirZ8);

	__m256 tminGTtzmax = _mm256_cmp_ps(tmin, tzmax, _CMP_GT_OS);
	__m256 tzminGTtmax = _mm256_cmp_ps(tzmin, tmax, _CMP_GT_OS);

	int checkZ = _mm256_movemask_ps(_mm256_or_ps(tminGTtzmax, tzminGTtmax));
	if (checkZ == 0)
		return false;

	return true;

	//from jacco's slides
	//__m128 t1 = _mm_mul_ps(_mm_sub_ps(bounds.bmin4, O4), rD4);
	//__m128 t2 = _mm_mul_ps(_mm_sub_ps(node->bmax4, O4), rD4);
	//__m128 vmax4 = _mm_max_ps(t1, t2), vmin4 = _mm_min_ps(t1, t2);
	//float* vmax = (float*)&vmax4, *vmin = (float*)&vmin4;
	//float tmax = min(vmax[0], min(vmax[1], vmax[2]));
	//float tmin = max(vmin[0], max(vmin[1], vmin[2]));
	//return tmax >= tmin && tmax >= 0;
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
