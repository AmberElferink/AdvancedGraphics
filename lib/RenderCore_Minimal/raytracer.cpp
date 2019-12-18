#include "core_settings.h"

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
bool Raytracer::IsOccluded( const Ray &ray, const Light &light )
{
	for ( const Mesh &mesh : scene.meshList )
	{
		int vertexCount = mesh.vcount / 3;
		for ( int i = 0; i < vertexCount; i++ )
		{
			Intersection intersection;
			if ( BVHNode::Intersect( ray, mesh.triangles[i], scene.matList, intersection) ) //If there are intersections
			{
				//light comes from infinetely far away
				if ( light.directionalLight )
					return true;
				else //light comes from a given point
				{
					if ( length( intersection.point - ray.O ) < length( light.position - ray.O ) ) //Between the light and the origin, not after
						return true;
				}
			}
		}
	}
	return false; //false if no intersections are found
}

/*Method that shoots a shadow ray and checks whether there are objects between the current intersection point and a light source*/
bool Raytracer::viewLight( Intersection intersection, const Light &light, float3 &lightVector )
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
bool Raytracer::viewDirLight( Intersection intersection, const Light &light, float3 &lightVector )
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
int Raytracer::viewSpotLight( Intersection intersection, const Light &light, float3 &lightVector )
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

/*method that checks whether a random point in an area of light is visible*/
bool Raytracer::viewAreaLight( const Intersection intersection, Light &light )
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

Intersection Raytracer::nearestIntersection( Ray ray )
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


float Raytracer::Fresnel(const float cosi, const float ncalc, const float n1, const float n2)
{
	//Fresnels law (how much reflects vs transmits).
//slide 20
//precalculations:
	float m = ncalc * sin(acos(cosi)); // (n1/n2) * sin(theta i)
	float coso = sqrtf(1 - (m*m));
	float n1i = n1 * cosi;
	float n2t = n2 * coso;
	float n1t = n1 * coso;
	float n2i = n2 * cosi;
	float sPolarisedRoot = (n1i - n2t) / (n1i + n2t);
	float pPolarisedRoot = (n1t - n2i) / (n1t + n2i);

	//refracted light percentage
	return 0.5f * (sPolarisedRoot * sPolarisedRoot + pPolarisedRoot * pPolarisedRoot);
}

int maxReflectionDepth = 7;
int reflectionDepth = 0;
//n1 default is air refraction index
//eta in lighthouse is 1/n of that material
float3 Raytracer::calcDielectric(Ray ray, Intersection intersection, const Intersection prevIntersection, int reflectionDepth, float n1)
{
	//Snells law:
	//formula: Advanced Graphics slides lecture 2 - Whitted Style - slide 18
	//or For a full derivation, see http://www.flipcode.com/archives/reflection_transmission.pdf

	float cosi = dot( intersection.norm, make_float3( -ray.D.x, -ray.D.y, -ray.D.z ) );
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
	float k = 1 - ( ncalc * ncalc ) * ( 1 - ( cosi * cosi ) );

	if (k < 0) //total internal reflection
	{
		return Reflect(ray, intersection, reflectionDepth); //make_float3(1, 0, 0);// 
	}


	float3 T = ncalc * ray.D + intersection.norm * (ncalc * cosi - sqrtf(k));
	Ray transmissionRay( intersection.point + 2 * EPSILON * T, T);

	float Fr = Fresnel(cosi, ncalc, n1, n2);
	float Ft = 1 - Fr; //transmitted light percentage

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

	float3 transmissionColor = Ft * Trace(transmissionRay, intersection, ++reflectionDepth);
	float3 reflectionColor = Fr * Reflect(ray, intersection, reflectionDepth); //Fr * make_float3(0, 0, 1);
	return transmissionColor + reflectionColor;
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
float3 Raytracer::Trace(const Ray &ray, const Intersection prevIntersection, int reflectionDepth)
{
	Intersection intersection = nearestIntersection( ray );

	if ( intersection.t > 10e29 )
	{
		reflectionDepth = -1;
		return make_float3(0.3, 0.3, 0.3) * ray.I; //background color
	}

	if ( reflectionDepth < maxReflectionDepth )
	{
		//Case of (partially) reflective material
		if (intersection.material.metallic)
			return Reflect(ray, intersection, reflectionDepth);
		else if (intersection.material.dielectric)
			return calcDielectric(ray, intersection, prevIntersection, reflectionDepth);
	}

	reflectionDepth = -1;
	//completely diffuse or maximum reflection depth
	return DirectIllumination(intersection) * ray.I; // ray.I is the intensity that comes through glass if it has passed through
}



float3 Raytracer::DirectIllumination( Intersection intersection )
{
	float3 intersectionColor = make_float3( 0, 0, 0 );
	float3 materialColor = intersection.material.color;
	if ( intersection.material.texture )
		{
		uint color = 0;
		TextureColor( intersection, intersection.triangle, color );
		materialColor.x = (float)((color >> 16) & 0x0000ff)/255.f;
		materialColor.y = (float)(( color >> 8 ) & 0x0000ff) / 255.f;
		materialColor.z = (float)(color & 0x0000ff) / 255.f;
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
				if (dotPr > 0)
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


void Raytracer::rayTrace( Bitmap *screen, const ViewPyramid &view, const int targetTextureID )
{
	for ( int j = 0; j < screen->height; j++ )
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
void Raytracer::rayTraceLine(Bitmap *screen, const ViewPyramid &view, const int targetTextureID, const int lineNr)
{
	int j = lineNr;
	for ( int i = 0; i < screen->width; i++ )
	{
		//u and v are the vectors within the virtual screen scaled between 0 and 1, so u = px / screenwidth and y = py / screenwidth
		float u = (float)i / (float)screen->width;
		float v = (float)j / (float)screen->height;
		float3 P = view.p1 + u * ( view.p2 - view.p1 ) + v * ( view.p3 - view.p1 );

		float3 dir = P - view.pos;		// vector in the direction you want to shoot your ray
		float3 D = dir / length( dir ); //normalize it

		Ray ray = Ray( view.pos, D );

		float3 intersectionColor = Trace(ray, Intersection(), 0);

		screen->pixels[i + j * screen->width] = FloatToIntColor( intersectionColor );
	}
}

void Raytracer::rayTraceBlock(const ViewPyramid &view, Bitmap *screen, const int targetTextureID, int lineStart, int lineEnd)
{

	for (int i =lineStart; i < lineEnd; i++)
	{
		rayTraceLine(screen, view, targetTextureID, i);
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
float3 prevp1 = make_float3(1000000); //random thingy
int framecounter = 0;

Bitmap* Raytracer::rayTraceRandom(const ViewPyramid &view, const int targetTextureID, int &frameCounter)
{
	//1 ray per light source geeft noisy image, maar andere random numbers = different numbers, optellen bij de accumulator. 
	//Die is 2x zo bright, maar /2 geeft weer normaal antwoord. Dat blijf je doen.

	storeBVH();

	if (!(view.p1.x == prevp1.x && view.p1.y == prevp1.y && view.p1.z == prevp1.z))
	{
		buffer->Clear();
		framecounter = 1;
	}

	for (int i = 0; i < 100; i++) 
	{

		//u and v are the vectors within the virtual screen scaled between 0 and 1, so u = px / screenwidth and y = py / screenwidth
		float u = RandomFloat();
		float v = RandomFloat();
		float3 P = view.p1 + u * (view.p2 - view.p1) + v * (view.p3 - view.p1);

		float3 dir = P - view.pos;		// vector in the direction you want to shoot your ray
		float3 D = dir / length(dir); //normalize it

		Ray ray = Ray(view.pos, D);

		float3 intersectionColor = Trace(ray, Intersection(), 0);
		int index = (int)((u * buffer->height) + (v * buffer->height) * buffer->width);
		buffer->pixels[index] += FloatToIntColor(intersectionColor);
	}
	framecounter++;
	return buffer;
}

void Raytracer::storeBVH()
{
	bvh.resize(scene.meshList.size());
	int i = 0;
	for (Mesh &mesh : scene.meshList)
	{
		bvh[i].ConstructBVH( mesh.vertices, mesh.vcount, mesh.triangles );
		i++;
	}
}