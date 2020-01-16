#include "core_settings.h"
int rayNr = 0;
// adapted from M�ller�Trumbore intersection algorithm: https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
/*bool Raytracer::Intersect( const Ray &ray, const CoreTri &triangle, Intersection &intersection )
{
	//TODO: het kan zijn dat een aantal dingen al geprecalculate zijn in CoreTri. Kijk daarnaar voor versnelling
	float3 vertex0 = triangle.vertex0;
	float3 vertex1 = triangle.vertex1;
	float3 vertex2 = triangle.vertex2;
	float3 edge1, edge2, h, s, q;
	float a, f, u, v;
	edge1 = vertex1 - vertex0;
	edge2 = vertex2 - vertex0;
	h = cross( ray.D, edge2 );
	a = dot( edge1, h );
	if ( a > -EPSILON && a < EPSILON )
		return false; // This ray is parallel to this triangle.
	f = 1.0f / a;
	s = ray.O - vertex0;
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
		intersection.material = *scene.matList[triangle.material];
		intersection.triangle = triangle;
		return true;
	}
	else // This means that there is a line intersection but not a ray intersection.
		return false;
}*/

/* Method that checks whether there are any objects between a light point and the origin of a shadowray */
//bool Raytracer::IsOccluded( const Ray &ray, const Light &light )
//{
//	for ( const Mesh &mesh : scene.meshList )
//	{
//		int vertexCount = mesh.vcount / 3;
//		for ( int i = 0; i < vertexCount; i++ )
//		{
//			Intersection intersection;
//			if ( BVHNode::Intersect( ray, mesh.triangles[i], scene.matList, intersection) ) //If there are intersections
//			{
//				//light comes from infinetely far away
//				if ( light.directionalLight )
//					return true;
//				else //light comes from a given point
//				{
//					if ( length( intersection.point - ray.O ) < length( light.position - ray.O ) ) //Between the light and the origin, not after
//						return true;
//				}
//			}
//		}
//	}
//	return false; //false if no intersections are found
//}

bool Raytracer::IsOccluded( const Ray &ray, const Light &light )
{
	int id = 0;
	Intersection intersection;
	for ( Mesh &mesh : scene.meshList )
	{
		if ( bvh[id].root->TraverseToFirst( ray, bvh[id].pool, bvh[id].indices, bvh[id].triangles, intersection, scene.matList ) ) // if there is at least one intersection
		{
			//light comes from infinitely far away
			if ( light.directionalLight )
				return true;
			else //light comes from a given point
			{
				if ( length( intersection.point - ray.O ) < length( light.position - ray.O ) ) //Between the light and the origin, not after
					return true;
			}
		}
		id++;
	}
	return false; //false if no intersections are found
}

/*Method that shoots a shadow ray and checks whether there are objects between the current intersection point and a light source*/
bool Raytracer::viewLight( const Intersection &intersection, const Light &light, float3 &lightVector )
{
	float3 dir = light.position - intersection.point; //vector between light and intersection point
	float dist = length( dir );
	lightVector = dir / dist; //normalized vector

	Ray shadowRay = Ray( intersection.point + intersection.norm * 0.0002f, lightVector ); //shadow ray from origin to light point

	if ( IsOccluded( shadowRay, light ) )
		return false; //cannot see light source
	else
		return true; //no objects that obstruct view of light source
}

/*Method that checks whether a directional light source can be viewed*/
bool Raytracer::viewDirLight( const Intersection &intersection, const Light &light, float3 &lightVector )
{
	lightVector = light.direction / length( light.direction ); //normalized vector

	Ray shadowRay = Ray( intersection.point + intersection.norm * 0.0002f, lightVector ); //shadow ray from origin to light point

	if ( IsOccluded( shadowRay, light ) )
		return false; //cannot see light source
	else
		return true; //no objects that obstruct view of light source
}

/*Method that checks whether the current intersection point can be seen from a spot light.
  It returns 0 if the spot light is occluded,
  1 if it can see the outer circle of light
  and 2 if it can also see the inner circle. */
int Raytracer::viewSpotLight( const Intersection &intersection, const Light &light, float3 &lightVector )
{
	//Normalized spotlight vector
	float3 spotDir = light.direction / length( light.direction );

	//Normalized vector from origin of the spot light to the intersection point
	float3 intersectionDir = light.position - intersection.point;
	intersectionDir = intersectionDir / length( intersectionDir );

	//Angle between the two vectors
	float cosInt = dot( spotDir, intersectionDir );

	//Point is not in the circle of light, so no need to check whether the point is occluded
	if ( cosInt < light.cosOuter )
		return false;

	lightVector = 1 * intersectionDir;

	Ray shadowRay = Ray( intersection.point + intersection.norm * 0.0002f, lightVector );

	if ( IsOccluded( shadowRay, light ) ) //not visible
		return 0;
	else if ( cosInt < light.cosInner ) //visible, but in outer circle
		return 1;
	else //in inner circle
		return 2;
}

/*Method that returns a random point on a triangle */
float3 Raytracer::randomPointTri( const CoreLightTri &triangle )
{
	//Generate two random floats between 0 and 1
	float u1 = ( (float)rand() ) / (float)RAND_MAX;
	float u2 = ( (float)rand() ) / (float)RAND_MAX;
	//Compute a random point on the quadrileteral that consists of the original triangle and the reflection of that triangle
	float3 point = triangle.vertex0 + u1 * ( triangle.vertex1 - triangle.vertex0 ) + u2 * ( triangle.vertex2 - triangle.vertex0 );

	//Compute a point p on the edge between vertex1 and vertex2 that forms a perpendicular line from the random point to p
	float3 A = point - triangle.vertex1;
	float3 normA = A / length( A );
	float3 edge = triangle.vertex2 - triangle.vertex1;
	float3 normEdge = edge / length( edge );

	float cosTheta = dot( normA, normEdge );

	float3 p = cosTheta * A;

	//compute the angle between the vector from p to the point and from p to vertex0
	float3 dirV0 = triangle.vertex0 - p;
	float3 normV0 = dirV0 / length( dirV0 );
	float3 dirA = A - p;
	float3 normVA = dirA / length( dirA );
	cosTheta = dot( normV0, normVA );

	//if this angle is larger than 90 degrees (and thus cosTheta < 0), then the random point is in the reflected triangle
	if ( cosTheta >= 0 )
		return point;
	else
		return point - 2 * p; //therefore, we reflect the point over p
}

/*Method that returns a random point on a random light*/
float3 Raytracer::randomPointLight( float3 &pl, Light &light_out )
{
	float rp = scene.areaLights * ( (float)rand() ) / (float)RAND_MAX;
	float area = 0;
	for (Light &light : scene.lightList)
	{
		if ( light.triangle.area + area >= rp )
			{
				light_out = light;
				pl = randomPointTri( light.triangle );
				return pl;
			}
			area += light.triangle.area;
	}
	return pl;
}

	/*method that checks whether a random point in an area of light is visible*/
	bool Raytracer::viewAreaLight( const Intersection &intersection, Light &light )
{
	//Random point on the triangle
	float3 point = randomPointTri( light.triangle );

	//normalized vector from intersection point to the random point on the triangle
	float3 dir = point - intersection.point;
	dir = dir / length( dir );

	Ray shadowRay = Ray( intersection.point + intersection.norm * 0.0002f, dir );

	//update position of current point on the light source
	light.position = point;

	return !IsOccluded( shadowRay, light ); //if not obstructed return true
}

uint Raytracer::FloatToIntColor( float3 floatColor )
{
	float r = min( floatColor.x, 1.0f );
	float g = min( floatColor.y, 1.0f );
	float b = min( floatColor.z, 1.0f );
	return ( ( ( uint )( r * 255.0f ) << 16 ) + ( ( uint )( g * 255.0f ) << 8 ) + ( uint )( b * 255.0f ) );
}

Intersection Raytracer::nearestIntersection( const Ray &ray )
{
	Intersection closest; //this will be your closest intersection of which you want to know the color
	closest.t = 10e30;
	closest.material = Material( make_float3( 0, 0, 0 ) ); //default black (background)
	int id = 0;
	Intersection save;
	save.t = 10e30;
	save.material = Material( make_float3( 0, 0, 0 ) ); //default black (background)

	//	Find closest intersection point for all meshes
	for ( Mesh &mesh : scene.meshList )
	{
		bvh[id].root->Traverse( ray, bvh[id].pool, bvh[id].indices, bvh[id].triangles, closest, scene.matList );
		if ( closest.t < save.t )
			save = closest;
		id++;
	}

	return save;
}

void Raytracer::nearestIntersection( const Ray8 &ray, Intersection8 &closest )
{
	int id = 0;

	//	Find closest intersection point for all meshes
	for ( Mesh &mesh : scene.meshList )
	{
		bvh[id].root->Traverse( ray, bvh[id].pool, bvh[id].indices, bvh[id].triangles, closest, scene.matList );
		id++;
	}
}

void Raytracer::nearestIntersection( Rays &r, Intersections &closests )
{
	int id = 0;

	//	Find closest intersection point for all meshes
	for ( Mesh &mesh : scene.meshList )
	{
		bvh[id].root->Traverse( r, RAYPACKETSIZE, Indices(), bvh[id].pool, bvh[id].indices, bvh[id].triangles, closests, scene.matList );
		id++;
	}
}

float Raytracer::Fresnel( const float cosi, const float ncalc, const float n1, const float n2 )
{
	//Fresnels law (how much reflects vs transmits).
	//slide 20
	//precalculations:
	float m = ncalc * sin( acos( cosi ) ); // (n1/n2) * sin(theta i)
	float coso = sqrtf( 1 - ( m * m ) );
	float n1i = n1 * cosi;
	float n2t = n2 * coso;
	float n1t = n1 * coso;
	float n2i = n2 * cosi;
	float sPolarisedRoot = ( n1i - n2t ) / ( n1i + n2t );
	float pPolarisedRoot = ( n1t - n2i ) / ( n1t + n2i );

	//refracted light percentage
	return 0.5f * ( sPolarisedRoot * sPolarisedRoot + pPolarisedRoot * pPolarisedRoot );
}

int maxReflectionDepth = 7;
int reflectionDepth = 0;
//n1 default is air refraction index
//eta in lighthouse is 1/n of that material
float3 Raytracer::calcDielectric( Ray ray, Intersection intersection, const Intersection prevIntersection, int reflectionDepth, float n1 )
{
	//Snells law:
	//formula: Advanced Graphics slides lecture 2 - Whitted Style - slide 18
	//or For a full derivation, see http://www.flipcode.com/archives/reflection_transmission.pdf

	float cosi = dot( intersection.norm, make_float3( -ray.D.x, -ray.D.y, -ray.D.z ) );
	float n2;
	if ( cosi <= 0 ) //you're going from n2 into n1, which makes them switch
	{
		n2 = n1;
		n1 = intersection.material.indexOfRefraction;
		cosi = -cosi;
		intersection.norm = -intersection.norm;
	}
	else //you're going from n1 into n2.
	{
		n2 = intersection.material.indexOfRefraction;
	}

	float ncalc = n1 / n2;

	//number within the root
	float k = 1 - ( ncalc * ncalc ) * ( 1 - ( cosi * cosi ) );

	if ( k < 0 ) //total internal reflection
	{
		return Reflect( ray, intersection, reflectionDepth ); //make_float3(1, 0, 0);//
	}

	float3 T = ncalc * ray.D + intersection.norm * ( ncalc * cosi - sqrtf( k ) );
	Ray transmissionRay( intersection.point + 2 * EPSILON * T, T );

	float Fr = Fresnel( cosi, ncalc, n1, n2 );
	float Ft = 1 - Fr; //transmitted light percentage

	//Beers law:
	if ( prevIntersection.material.dielectric ) //there was a previous intersection in dielectric
	{
		float distance = length( intersection.point - prevIntersection.point );
		float3 I;
		I.x = exp( -intersection.material.absorption.x * distance );
		I.y = exp( -intersection.material.absorption.y * distance );
		I.z = exp( -intersection.material.absorption.z * distance );
		transmissionRay.I *= I;
		ray.I *= I;
	}

	float3 transmissionColor = Ft * Trace( transmissionRay, intersection, ++reflectionDepth );
	float3 reflectionColor = Fr * Reflect( ray, intersection, reflectionDepth ); //Fr * make_float3(0, 0, 1);
	return transmissionColor + reflectionColor;
}

float3 Raytracer::DielectricPath(Ray ray, Intersection intersection, const Intersection prevIntersection, float n1)
{
	//Snells law:
	//formula: Advanced Graphics slides lecture 2 - Whitted Style - slide 18
	//or For a full derivation, see http://www.flipcode.com/archives/reflection_transmission.pdf

	float cosi = dot(intersection.norm, make_float3(-ray.D.x, -ray.D.y, -ray.D.z));
	float n2;
	if (cosi <= 0) 	 //you're going from n2 into n1, which makes them switch	
	{
		n2 = n1;
		n1 = intersection.material.indexOfRefraction;
		cosi = -cosi;
		intersection.norm = -intersection.norm;
	}
	else //you're going from n1 into n2.
	{
		n2 = intersection.material.indexOfRefraction;
	}

	float ncalc = n1 / n2;

	//number within the root
	float k = 1 - (ncalc * ncalc) * (1 - (cosi * cosi));

	if (k < 0) //total internal reflection
	{
		return Sample(ray.Reflect(intersection.norm, intersection.point), intersection, reflectionDepth); //make_float3(1, 0, 0);// 
	}


	float3 T = ncalc * ray.D + intersection.norm * (ncalc * cosi - sqrtf(k));
	Ray transmissionRay(intersection.point + 2 * EPSILON * T, T);

	float rn = (float)rand() / (float)RAND_MAX;
	float Fr = Fresnel(cosi, ncalc, n1, n2);
	if (rn < Fr) //reflect
		return Sample(ray.Reflect(intersection.norm, intersection.point), intersection, reflectionDepth);
	else //transmit and apply beers law
	{
		//Beers law:
		if (prevIntersection.material.dielectric) //there was a previous intersection in dielectric
		{
			float distance = length(intersection.point - prevIntersection.point);
			float3 I;
			I.x = exp(-intersection.material.absorption.x * distance);
			I.y = exp(-intersection.material.absorption.y * distance);
			I.z = exp(-intersection.material.absorption.z * distance);
			transmissionRay.I *= I;
			ray.I *= I;
		}
		return Sample(transmissionRay, intersection, ++reflectionDepth);
	}

}

void Raytracer::TextureColor( Intersection &intersection, const CoreTri &triangle, uint &color )
{
	Texture tex = Texture( intersection.material.texture->width, intersection.material.texture->height );

	int height = tex.height;
	int width = tex.width;

	//normalized vector from v0 to p
	float3 S = intersection.point - triangle.vertex0;
	float3 normS = S / length( S );

	//normalized vector from v0 to the point closest to p on edge(v0,v1)
	float3 U = triangle.vertex1 - triangle.vertex0;
	float3 normU = U / length( U );
	float lengthU = length( U );

	//normalized vector from v0 to the point closest to p on edge(v0,v2)
	float3 V = triangle.vertex2 - triangle.vertex0;
	float3 normV = V / length( V );
	float lengthV = length( V );

	//respective angles with vector (v0,p)
	float cosX = dot( normS, normU );
	float cosY = dot( normS, normV );

	float u = length( intersection.point - triangle.vertex0 ) * cosX;
	float v = length( intersection.point - triangle.vertex0 ) * cosY;

	uint totalPixels = width * height;
	float scaleX = (float)width / length( triangle.vertex1 - triangle.vertex0 );
	float scaleY = (float)height / length( triangle.vertex2 - triangle.vertex0 );

	//if ( intersection.point.x < triangle.v0 && intersection.point.y < triangle.v0 )
	//{
	//	u = lengthU - u;
	//	v = lengthV - v;
	//}

	uint index = ( v * scaleY * width + u * scaleX );

	color = intersection.material.texture->pixels[index];
}

float3 Raytracer::Reflect( const Ray &ray, const Intersection &intersection, int reflectionDepth )
{
	//s denotes the amount of light that is reflected and d the amount that is absorbed
	float s = intersection.material.specularity;

	if ( s == 0 ) //no reflection
		return DirectIllumination( intersection );

	float d = 1 - intersection.material.specularity;

	//Computes the direction of the reflected ray
	float3 reflectedDir = ray.D - 2 * dot( ray.D, intersection.norm ) * intersection.norm;
	Ray reflectedRay = Ray( intersection.point + 2 * EPSILON * reflectedDir, reflectedDir );

	if ( d == 0 ) //no absorption
		return Trace( reflectedRay, intersection, ++reflectionDepth );
	else
	{
		return s * intersection.material.color * Trace( reflectedRay, intersection, ++reflectionDepth ) + d * DirectIllumination( intersection );
	}
}



//Method that sends a ray into a scene and returns the color of the hitted objects
//prevIntersection is only used for dieelectric n2.
Color8 Raytracer::Trace( Ray8 &ray, const Intersection8 prevIntersection, int reflectionDepth )
{
	Intersection8 closest;
	nearestIntersection( ray, closest );

	for ( int i = 0; i < 8; i++ )
	{
		if ( closest.t[i] > 10e29 ) //background
		{
			reflectionDepth = -1;
			if ( isnan( abs( ray.deadMask[i] ) ) ) //ray.activeMask[i] > 0)
			{
				ray.color.b[i] = 0.3f; //background color //should be * ray.I
				ray.color.g[i] = 0.3f;
				ray.color.r[i] = 0.3f;
				ray.deadMask[i] = 0x00000000;
			}
		}
		else //direct illumination
		{
			if ( isnan( abs( ray.deadMask[i] ) ) ) //it's -nan, which means it's true. (0 is false)
			{
				float3 c1 = DirectIllumination( closest.intersections[i] ); //should be * ray.I
				ray.color.b[i] = c1.x;
				ray.color.g[i] = c1.y;
				ray.color.r[i] = c1.z;
				ray.deadMask[i] = 0x00000000;
			}
		}
		if ( _mm256_movemask_ps( ray.activeMask8 ) == 0 )
		{
			return ray.color;
		}
	}
}

//Method that sends a ray into a scene and returns the color of the hitted objects
//prevIntersection is only used for dieelectric n2.
void Raytracer::Trace( Rays &r, Indices I, const Intersections prevIntersection, int reflectionDepth )
{
	Intersections closests;
	nearestIntersection( r, closests );

	//if (reflectionDepth == 0) //color the background for the first round for all inactive rays
	//{
	//	for (int j = r.ia; j < RAYPACKETSIZE; j++)
	//	{
	//		for (int i = 0; i < 8; i++)
	//		{
	//			//inactive rays got the background
	//			r.rays[j].color.r[i] = 0.3f; //background color //should be * ray.I
	//			r.rays[j].color.g[i] = 0.3f;
	//			r.rays[j].color.b[i] = 0.3f;
	//			r.rays[j].deadMask[i] = 0x00000000;
	//		}
	//	}
	//}

	//if (r.ia == 0)
	//		return; //there are no rays that still need a color

	//all active rays must be checked.
	for ( int j = 0; j < RAYPACKETSIZE; j++ )
	{
		for ( int i = 0; i < 8; i++ )
		{
			if ( closests.inter[I.I[j]].t[i] < 10e29 ) //direct illumination
			{
				if ( isnan( abs( r.rays[I.I[j]].deadMask[i] ) ) ) //it's -nan, which means it's true. (0 is false)
				{
					float3 c1 = DirectIllumination( closests.inter[I.I[j]].intersections[i] ); //should be * ray.I
					r.rays[I.I[j]].color.b[i] = c1.x;
					r.rays[I.I[j]].color.g[i] = c1.y;
					r.rays[I.I[j]].color.r[i] = c1.z;
					r.rays[I.I[j]].deadMask[i] = 0x00000000;
				}
			}
			else //background
			{
				reflectionDepth = -1;
				if ( isnan( abs( r.rays[j].deadMask[i] ) ) ) //ray.activeMask[i] > 0)
				{
					r.rays[I.I[j]].color.b[i] = 0.3f; //background color //should be * ray.I
					r.rays[I.I[j]].color.g[i] = 0.3f;
					r.rays[I.I[j]].color.r[i] = 0.3f;
					r.rays[I.I[j]].deadMask[i] = 0x00000000;
				}
			}
			//if (_mm256_movemask_ps(r.rays[r.I[j]].activeMask8) == 0)
			//{
			//	swap(rI[j], r.I[r.ia]);
			//	r.ia--;
			//}
		}
	}

	//if ( reflectionDepth < maxReflectionDepth )
	//{
	//	//Case of (partially) reflective material
	//	if (intersection.material.metallic)
	//		return Reflect(ray, intersection, reflectionDepth);
	//	else if (intersection.material.dielectric)
	//		return calcDielectric(ray, intersection, prevIntersection, reflectionDepth);
	//}

	//reflectionDepth = -1;

	//completely diffuse or maximum reflection depth
	//return DirectIllumination(intersection) * ray.I; // ray.I is the intensity that comes through glass if it has passed through
}

//Method that sends a ray into a scene and returns the color of the hitted objects
//prevIntersection is only used for dieelectric n2.
float3 Raytracer::Trace( const Ray &ray, const Intersection prevIntersection, int reflectionDepth )
{
	Intersection intersection = nearestIntersection( ray );

	if ( intersection.t > 10e29 )
	{
		reflectionDepth = -1;
		return make_float3( 0.3, 0.3, 0.3 ) * ray.I; //background color
	}

	if ( reflectionDepth < maxReflectionDepth )
	{
		//Case of (partially) reflective material
		if ( intersection.material.metallic )
			return Reflect( ray, intersection, reflectionDepth );
		else if ( intersection.material.dielectric )
			return calcDielectric( ray, intersection, prevIntersection, reflectionDepth );
	}

	reflectionDepth = -1;
	//completely diffuse or maximum reflection depth
	float3 color = DirectIllumination( intersection ) * ray.I; // ray.I is the intensity that comes through glass if it has passed through
	return color;
}

float3 Raytracer::DirectIllumination( Intersection intersection )
{
	float3 intersectionColor = make_float3( 0, 0, 0 );
	float3 materialColor = intersection.material.color;
	if ( 0 ) //intersection.material.texture )
	{
		uint color = 0;
		TextureColor( intersection, intersection.triangle, color );
		materialColor.x = (float)( ( color >> 16 ) & 0x0000ff ) / 255.f;
		materialColor.y = (float)( ( color >> 8 ) & 0x0000ff ) / 255.f;
		materialColor.z = (float)( color & 0x0000ff ) / 255.f;
	}

	/*Check for all lights whether they can be seen from the current intersection point*/
	for ( Light &light : scene.lightList )
	{
		//std::cout << closest.material.diffuse.x << " " << closest.material.diffuse.y << " " << closest.material.diffuse.z << endl;
		float3 lightVector;
		if ( light.pointLight )
		{
			if ( viewLight( intersection, light, lightVector ) )
			{
				float dist = length( intersection.point - light.position );
				float dotPr = dot( intersection.norm, lightVector );
				if ( dotPr > 0 )
				{
					intersectionColor += materialColor * light.radiance * ( 1 / ( dist * dist ) ) * dotPr; //If light source can be seen, multiply color with current pixel color
					int w = 0;
				}
			}
		}
		else if ( light.directionalLight )
		{
			if ( viewDirLight( intersection, light, lightVector ) )
			{
				float dotPr = dot( intersection.norm, lightVector );
				if ( dotPr > 0 )
					intersectionColor += materialColor * light.radiance * dotPr;
			}
		}
		else if ( light.spotLight )
		{
			int option = viewSpotLight( intersection, light, lightVector );
			if ( option != 0 ) //spotlight is visible
			{
				float dist = length( light.position - intersection.point );
				float dotPr = dot( intersection.norm, lightVector );

				//set the brightness difference between the inner and outer circle of the spotlight
				float difference = 1;
				if ( option == 2 )
					difference += 0.5;
				if ( dotPr > 0 )
					intersectionColor += materialColor * light.radiance * dotPr * difference / ( dist * dist );
			}
		}
		else //(if light.areaLight)
		{
			//Send a number of random rays to the areaLight
			int visible = 0;
			int k = 20;
			for ( int i = 0; i < k; i++ )
			{
				if ( viewAreaLight( intersection, light ) )
					visible += 1;
			}
			float dist = length( light.triangle.centre - intersection.point );
			float3 lightVector = ( light.triangle.centre - intersection.point ) / dist;
			float dotPr = dot( intersection.norm, lightVector );

			if ( dotPr > 0 )
				intersectionColor += materialColor * light.triangle.radiance * light.triangle.area * dotPr * visible / (float)k / ( dist * dist );
		}
	}
	return intersectionColor;
}

//---------------------------------------------------------PATHTRACER METHODS------------------------------------

/* Random vector on the hemisphere following a uniform distribution */
float3 Raytracer::DiffuseReflection( float3 N )
{
	while ( true )
	{
		//Generate three random floats between -1 and 1
		float x = 2 * ( (float)rand() ) / (float)RAND_MAX;
		float y = 2 * ( (float)rand() ) / (float)RAND_MAX;
		float z = 2 * ( (float)rand() ) / (float)RAND_MAX;
		if ( x > 1 )
			x = 1 - x;
		if ( y > 1 )
			y = 1 - y;
		if ( z > 1 )
			z = 1 - z;

		float radius = x * x + y * y + z * z;
		if ( radius <= 1 )
		{
			float3 dir = make_float3( x, y, z );
			dir = dir / sqrt( radius );
			if ( dot( N, dir ) < 0 )
				dir = -dir;
			return dir;
		}
	}
}

/* Random vector on the hemisphere following a cosine weighted distribution */
float3 Raytracer::CosineWeightedDiffuseReflection( float3 N )
{
	float r0 = (float)rand() / (float)RAND_MAX;
	float r1 = (float)rand() / (float)RAND_MAX;
	float r = sqrt( r0 );
	float theta = 2 * PI * r1;
	float x = r * cosf( theta );
	float y = r * sinf( theta );
	float3 P = make_float3( x, y, sqrt( 1 - r0 ) );

	//Vector P needs to be transformed to tangent space
	float3 W;
	if ( abs( N.x ) > 0.99 )
		W = make_float3( 0, 1, 0 );
	else
		W = make_float3( 1, 0, 0 );
	float3 T = cross( N, W );
	T = T / length(T);
	float3 B = cross( T, N );
	//return make_float3( dot( P, T ), dot( P, B ), dot( P, N ) );
	return P.x * T + P.y * B + P.z * N;
}
float3 Raytracer::MISample(Ray &ray, bool lastSpecular)
{
	float3 T = make_float3(1);
	float3 E = make_float3(0);
	while ( true )
	{
		Intersection I = nearestIntersection( ray );
		float3 BRDF = I.material.color / PI;
		if ( I.t > 10e29 ) 
			return E;
		if ( I.triangle.ltriIdx >= 0 )
			{
			if ( lastSpecular )
				return I.material.color;
			return E;
			}
		//specular surfaces
		if ( I.material.metallic )
		{
			//take a random number between 0 and 1
			float rn = (float)rand() / (float)RAND_MAX;
			if ( rn <= I.material.specularity )
			{
				float3 reflectedDir = ray.D - 2 * dot( ray.D, I.norm ) * I.norm;
				Ray reflectedRay = Ray( I.point + 2 * EPSILON * reflectedDir, reflectedDir );
				return E + T * MISample( reflectedRay, true ) * I.material.color;
			}
		}
		// sample a random light source
		float3 pl; 
		Light light;
		randomPointLight(pl,light);
		float dist = length( pl - I.point );
		pl = ( pl - I.point ) / dist;
		Ray lr( I.point, pl );
		float dot1 = dot( I.norm, pl );
		float dot2 = dot( light.triangle.N, -pl );
		if ( dot1 > 0 && dot2 > 0 )
			{
				if ( !IsOccluded( lr, light ) )
				{
					float solidAngle = ( dot2 * light.triangle.area ) / ( dist * dist );
					float misPDF = 1 / solidAngle + 1 / (2*PI); //1 /(2*PI);
					E += T * ( dot1 / misPDF ) * light.triangle.radiance * BRDF;
				}
			}
		// continue random walk
		float3 R = CosineWeightedDiffuseReflection( I.norm );
		float dotR = dot( R, I.norm );
		float PDF = dotR / PI;
		Ray r( I.point, R );
		T *= ( dotR / PDF ) * BRDF;
		ray = r;
		lastSpecular = false;
	}
	return E;
}

float3 Raytracer::Sample(const Ray &ray, Intersection prevIntersection, int reflectionDepth)
{
	Intersection I = nearestIntersection( ray );
	// terminate if ray left the scene
	if ( I.t > 10e29 )
		return make_float3( 0 );

	//we do not want to sample light sources on our random walk
	if ( I.triangle.ltriIdx >= 0 )
		{
		if ( prevIntersection.material.metallic )
			return I.material.color;
		else return make_float3(0);
		}
	//specular surfaces
	if ( I.material.metallic )
	{
		//take a random number between 0 and 1
		float rn = (float)rand() / (float)RAND_MAX;
		if ( rn <= I.material.specularity )
		{
			return Sample( ray.Reflect( I.norm, I.point), I, reflectionDepth ) * I.material.color * ray.I;
		}
	}
	if (I.material.dielectric)
	{
		return DielectricPath(ray, I, prevIntersection);
	}

	float3 pl; //random point on random light
	float3 Ld = make_float3( 0 );
	float3 BRDF = I.material.color / PI;
	
	Light light;
	randomPointLight( pl, light );

	//normalize vector to the random point on a random light
	float dist = length( pl - I.point );
	pl = ( pl - I.point ) / dist;
	Ray lr( I.point, pl ); //send ray to random point on one of the lights

	float dot1 = dot( I.norm, pl );
	float dot2 = dot( light.triangle.N, -pl );
	if ( dot1 > 0 && dot2 > 0)
	{
		if (!IsOccluded(lr,light))
		{
			float solidAngle = ( dot2 * light.triangle.area ) / ( dist * dist );
			//USE FOR Multiple Importance Sampling: 
			//float misPDF = 1 / solidAngle + 1 / ( 2 * PI ); //1 /(2*PI);
			//Ld = ( dot1 / misPDF ) * light.triangle.radiance * BRDF;
			Ld = solidAngle * dot1 * light.triangle.radiance * BRDF;
		}
	}

	// continue in random direction
	float3 R = CosineWeightedDiffuseReflection( I.norm );
	Ray r( I.point, R );
	// update throughput
	float dotR = dot( I.norm, R );
	float PDF =  dotR / PI; //set to 1 / (2 * PI) if you are using DiffReflection
	float3 Ei = Sample( r, I, reflectionDepth ) * dotR / PDF * ray.I;
	return BRDF * Ei + Ld;
}

void Raytracer::rayTrace( Bitmap *screen, const ViewPyramid &view, const int targetTextureID )
{
	for ( uint j = 0; j < screen->height; j++ )
	{
		rayTraceLine( screen, view, targetTextureID, j );
	}
}

//-----------------------------------------------------
// shoot a ray through point p to intersect the scene.
//  p3 |------------------|
//     |        .         |
//     |        p         | screenheight
//  p1 |------------------| p2
//   screenwidth
//-----------------------------------------------------
void Raytracer::rayTraceLine( Bitmap *screen, const ViewPyramid &view, const int targetTextureID, const int lineNr )
{
	int j = lineNr;
	for ( uint i = 0; i < screen->width; i++ )
	{
		//printf("pX: %i, i: %i, pY: %i, j: %i\n", probePos.x, i, probePos.y, j);
		if ( i == probePos.x && j == probePos.y )
			printf( "probing: x: %i, y: %i, rayNr: %i\n", i, j, rayNr );

		//u and v are the vectors within the virtual screen scaled between 0 and 1, so u = px / screenwidth and y = py / screenwidth
		float u = (float)i / (float)screen->width;
		float v = (float)j / (float)screen->height;
		float3 P = view.p1 + u * ( view.p2 - view.p1 ) + v * ( view.p3 - view.p1 );

		float3 dir = P - view.pos;		// vector in the direction you want to shoot your ray
		float3 D = dir / length( dir ); //normalize it

		Ray ray = Ray( view.pos, D );

		float3 intersectionColor = Trace( ray, Intersection(), 0 );

		screen->pixels[i + j * screen->width] = FloatToIntColor( intersectionColor );
		rayNr++;
	}
}

void Raytracer::pathTrace( Bitmap *screen, const ViewPyramid &view, const int targetTextureID, uint sampleCount )
{
	for ( uint j = 0; j < screen->height; j++ )
	{
		for ( uint i = 0; i < screen->width; i++ )
		{
			//printf("pX: %i, i: %i, pY: %i, j: %i\n", probePos.x, i, probePos.y, j);
			if ( i == probePos.x && j == probePos.y )
				printf( "probing: x: %i, y: %i, rayNr: %i\n", i, j, rayNr );

			//u and v are the vectors within the virtual screen scaled between 0 and 1, so u = px / screenwidth and y = py / screenwidth
			float u = (float)i / (float)screen->width;
			float v = (float)j / (float)screen->height;
			float3 P = view.p1 + u * ( view.p2 - view.p1 ) + v * ( view.p3 - view.p1 );

			float3 dir = P - view.pos;		// vector in the direction you want to shoot your ray
			float3 D = dir / length( dir ); //normalize it

			Ray ray = Ray( view.pos, D );

			float3 intersectionColor = Sample(ray, Intersection(), 0);
			buffer->pixels[i + j * buffer->width] = (buffer->pixels[i + j * buffer->width] * (float) (sampleCount) + intersectionColor) / (float) (sampleCount + 1);
			screen->pixels[i + j * screen->width] = FloatToIntColor(buffer->pixels[i + j * buffer->width]);
			rayNr++;
		}
	}
	//rayNr = 0;
}

void Raytracer::rayTraceLineAVX( Bitmap *screen, const ViewPyramid &view, const int targetTextureID, const int lineNr )
{
	for ( uint i = 0; i < screen->width; i += 8 )
	{
		float3 O[8], D[8];

		for ( int di = 0; di < 8; di++ )
		{
			int x = i + di;
			//u and v are the vectors within the virtual screen scaled between 0 and 1, so u = px / screenwidth and y = py / screenwidth
			float u = ( (float)( i + di ) ) / (float)screen->width;
			float v = ( (float)lineNr ) / (float)screen->height;
			float3 P = view.p1 + u * ( view.p2 - view.p1 ) + v * ( view.p3 - view.p1 );

			float3 dir = P - view.pos; // vector in the direction you want to shoot your ray
			O[di] = view.pos;
			D[di] = dir / length( dir ); //normalize it
		}

		Color8 intersectionColor = Trace( Ray8( O, D ), Intersection8(), 0 );
		for ( int di = 0; di < 8; di++ )
		{
			screen->pixels[i + di + lineNr * screen->width] = FloatToIntColor( make_float3( intersectionColor.b[di], intersectionColor.g[di], intersectionColor.r[di] ) );
		}
		rayNr += 8;
	}
}

void Raytracer::rayTraceInPackets( Bitmap *screen, const ViewPyramid &view, const int targetTextureID, const int lineNr )
{
	for ( uint i = 0; i < screen->width; i += 8 )
	{
		Rays r;
		float3 O[8], D[8];
		for ( int dj = 0; dj < RAYPACKETSIZE; dj++ )
		{
			for ( int di = 0; di < 8; di++ )
			{
				if ( i + di == probePos.x && lineNr + dj == probePos.y )
					printf( "probing: x: %i, y: %i, rayNr: %i\n", i + di, lineNr + dj, rayNr );
				int x = i + di;
				//u and v are the vectors within the virtual screen scaled between 0 and 1, so u = px / screenwidth and y = py / screenwidth
				float u = ( (float)( i + di ) ) / (float)screen->width;
				float v = ( (float)( lineNr + dj ) ) / (float)screen->height;
				float3 P = view.p1 + u * ( view.p2 - view.p1 ) + v * ( view.p3 - view.p1 );

				float3 dir = P - view.pos; // vector in the direction you want to shoot your ray
				O[di] = view.pos;
				D[di] = dir / length( dir ); //normalize it
			}
			r.rays[dj] = Ray8( O, D );
		}

		Trace( r, Indices(), Intersections(), 0 );

		for ( int dj = 0; dj < RAYPACKETSIZE; dj++ )
		{
			for ( int di = 0; di < 8; di++ )
			{
				screen->pixels[( i + di ) + ( lineNr + dj ) * screen->width] = FloatToIntColor( make_float3( r.rays[dj].color.b[di], r.rays[dj].color.g[di], r.rays[dj].color.r[di] ) );
			}
		}
		rayNr += 8; //*RAYPACKETSIZE, but for debugging only look at the first row.
	}
}

void Raytracer::rayTraceBlock( const ViewPyramid &view, Bitmap *screen, const int targetTextureID, int lineStart, int lineEnd )
{

	for ( int i = lineStart; i < lineEnd; i++ )
	{
		rayTraceLine( screen, view, targetTextureID, i );
	}
}

void Raytracer::rayTraceBlockAVX( const ViewPyramid &view, Bitmap *screen, const int targetTextureID, int lineStart, int lineEnd )
{

	for ( int i = lineStart; i < lineEnd; i++ )
	{
		rayTraceLineAVX( screen, view, targetTextureID, i );
	}
}

//-----------------------------------------------------
// shoot a ray through point p to intersect the scene.
//  p3 |------------------|
//     |        .         |
//     |        p         | screenheight
//  p1 |------------------| p2
//   screenwidth
//-----------------------------------------------------
//float3 prevp1 = make_float3(1000000); //random thingy
//int framecounter = 0;
//
//Bitmap* Raytracer::rayTraceRandom(const ViewPyramid &view, const int targetTextureID, int &frameCounter)
//{
//	//1 ray per light source geeft noisy image, maar andere random numbers = different numbers, optellen bij de accumulator.
//	//Die is 2x zo bright, maar /2 geeft weer normaal antwoord. Dat blijf je doen.
//
//	storeBVH();
//
//	if (!(view.p1.x == prevp1.x && view.p1.y == prevp1.y && view.p1.z == prevp1.z))
//	{
//		buffer->Clear();
//		framecounter = 1;
//	}
//
//	for (int i = 0; i < 100; i++)
//	{
//
//		//u and v are the vectors within the virtual screen scaled between 0 and 1, so u = px / screenwidth and y = py / screenwidth
//		float u = RandomFloat();
//		float v = RandomFloat();
//		float3 P = view.p1 + u * (view.p2 - view.p1) + v * (view.p3 - view.p1);
//
//		float3 dir = P - view.pos;		// vector in the direction you want to shoot your ray
//		float3 D = dir / length(dir); //normalize it
//
//		Ray ray = Ray(view.pos, D);
//
//		float3 intersectionColor = Trace(ray, Intersection(), 0);
//		int index = (int)((u * buffer->height) + (v * buffer->height) * buffer->width);
//		buffer->pixels[index] += FloatToIntColor(intersectionColor);
//	}
//	framecounter++;
//	return buffer;
//}

void Raytracer::storeBVH()
{
	bvh.resize( scene.meshList.size() );
	int i = 0;
	for ( Mesh &mesh : scene.meshList )
	{
		bvh[i].ConstructBVH( mesh.vertices, mesh.vcount, mesh.triangles );
		i++;
	}
}