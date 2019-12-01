#include "core_settings.h"
#include <iostream>
// adapted from M�ller�Trumbore intersection algorithm: https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
bool Raytracer::Intersect(const Ray &ray, const CoreTri &triangle, Intersection &intersection)
{
	//TODO: het kan zijn dat een aantal dingen al geprecalculate zijn in CoreTri. Kijk daarnaar voor versnelling
	float3 vertex0 = triangle.vertex0;
	float3 vertex1 = triangle.vertex1;
	float3 vertex2 = triangle.vertex2;
	float3 edge1, edge2, h, s, q;
	float a, f, u, v;
	edge1 = vertex1 - vertex0;
	edge2 = vertex2 - vertex0;
	h = cross(ray.E, edge2);
	a = dot(edge1, h);
	if (a > -EPSILON && a < EPSILON)
		return false; // This ray is parallel to this triangle.
	f = 1.0 / a;
	s = ray.O - vertex0;
	u = f * dot(s, h);
	if (u < 0.0 || u > 1.0)
		return false;
	q = cross(s, edge1);
	v = f * dot(ray.E, q);
	if (v < 0.0 || u + v > 1.0)
		return false;
	// At this stage we can compute t to find out where the intersection point is on the line.
	float t = f * dot(edge2, q);
	if (t > EPSILON && t < 1 / EPSILON) // ray intersection
	{
		float3 intersectionPoint = ray.O + ray.E * t;
		float3 normal = make_float3(triangle.Nx, triangle.Ny, triangle.Nz); //TODO maybe use mesh N
		intersection = Intersection(t, intersectionPoint, normal, triangle);
		intersection.material = *scene.matList[triangle.material];
		return true;
	}
	else // This means that there is a line intersection but not a ray intersection.
		return false;
}

/* Method that checks whether there are any objects between a light point and the origin of a shadowray */
bool Raytracer::IsOccluded( const Ray &ray, const Light &light)
{
	for ( Mesh &mesh : scene.meshList )
	{
		int vertexCount = mesh.vcount / 3;
		for ( int i = 0; i < vertexCount; i++ )
		{
			Intersection intersection;
			if ( Intersect( ray, mesh.triangles[i], intersection ) ) //If there are intersections
				if (length(intersection.point - ray.O) <  length(light.position - ray.O)) //Between the light and the origin, not after
					return true;
		}
	}
	return false; //false if no intersections are found
}

/*Method that shoots a shadow ray and checks whether there are objects between the current intersection point and a light source*/
bool Raytracer::viewLight( Intersection intersection, const Light &light, float3 &lightVector )
{
	float3 dir = intersection.point - light.position; //vector between light and intersection point
	float dist = length(dir);
	lightVector = dir / dist; //normalized vector


	Ray shadowRay = Ray( intersection.point + lightVector * 0.0002f, lightVector ); //shadow ray from origin to light point

	Intersection closest;
	
	if ( IsOccluded( shadowRay, light ) ) 
		return false; //cannot see light source
	else
		return true; //no objects that obstruct view of light source
}

uint Raytracer::FloatToIntColor( float3 floatColor )
{
	float r = min(floatColor.x, 1.0f);
	float g = min(floatColor.y, 1.0f);
	float b = min(floatColor.z, 1.0f);
	return ( ((uint)( r * 255.0f ) << 16) + ((uint)(g * 255.0f ) << 8) + (uint)( b * 255.0f ) );
	}

Intersection Raytracer::nearestIntersection(Ray ray)
{
	Intersection closest; //this will be your closest intersection of which you want to know the color
	closest.t = 10e30;
	closest.material = Material( make_float3( 0, 0, 0 ) ); //default black (background)

	//Find closest intersection point for all meshes
	for ( Mesh &mesh : scene.meshList )
	{
		int triangleCount = mesh.vcount / 3;
		for ( int i = 0; i < triangleCount; i++ ) //find the closest triangle intersection for all triangles
		{
			Intersection intersection;
			if ( Intersect( ray, mesh.triangles[i], intersection ) )
			{
				if ( intersection.t < closest.t ) //update closest intersection point only if distance to intersection point is shorter
					closest = intersection;
			}
		}
	}

	return closest;
}

int maxReflectionDepth = 10;
int reflectionDepth = -1; //start at -1, the first trace is no reflection
//Method that sends a ray into a scene and returns the color of the hitted objects
float3 Raytracer::Trace(Ray ray)
{
	reflectionDepth++;

	Intersection intersection = nearestIntersection( ray );

	if (intersection.t > 10e29)
	{
		reflectionDepth = 0;
		return make_float3(0.2, 0.2, 0.2); //background color
	}

	//Case of (partially) reflective material
	if (intersection.material.metallic)
	{
		//s denotes the amount of light that is reflected and d the amount that is absorbed
		float s = intersection.material.specularity;

		if (s == 0) //no reflection
		{
			reflectionDepth = -1;
			return TotalLight( intersection );
		}

		float d = 1 - intersection.material.specularity;

		//Computes the direction of the reflected ray
		float3 reflectedDir = ray.E - 2 * dot(ray.E, intersection.norm) * intersection.norm;
		Ray reflectedRay = Ray(intersection.point + EPSILON * reflectedDir, reflectedDir);

		if (reflectionDepth < maxReflectionDepth)
		{
			if (d == 0) //no absorption
				return Trace(reflectedRay);
			else
				return s * intersection.material.diffuse * Trace(reflectedRay) + d * TotalLight(intersection);
		}
	}
	



	reflectionDepth = -1;
	//completely diffuse or maximum reflection depth
	return TotalLight(intersection);
}

float3 Raytracer::TotalLight(Intersection intersection)
{
	float3 intersectionColor = make_float3( 0, 0, 0 );

	/*Check for all lights whether they can be seen from the current intersection point*/
	for ( Light &light : scene.lightList )
	{
		//std::cout << closest.material.diffuse.x << " " << closest.material.diffuse.y << " " << closest.material.diffuse.z << endl;
		float3 lightVector;
		if ( viewLight( intersection, light, lightVector ) )
		{
			float dist = length( light.position - intersection.point );
			float dotPr = dot(intersection.norm, lightVector);
			if(dotPr > 0)
				intersectionColor += intersection.material.diffuse * light.radiance * (1 / ( dist * dist )) * dotPr; //If light source can be seen, multiply color with current pixel color
		}
	}
	return intersectionColor;
}

//-----------------------------------------------------
// shoot a ray through point p to intersect the scene.
//  p3 |------------------| 
//     |        .         |
//     |        p         | screenheight
//  p1 |------------------| p2
//   screenwidth
//-----------------------------------------------------
void Raytracer::rayTrace(Bitmap *screen, const ViewPyramid &view, const int targetTextureID)
{
	Timer t;
	for ( int j = 0; j < screen->height; j++ )
	{
		rayTraceLine(screen, view, targetTextureID, j);
	}
}

void Raytracer::rayTraceLine(Bitmap *screen, const ViewPyramid &view, const int targetTextureID, const int lineNr)
{
	int j = lineNr;
	for (int i = 0; i < screen->width; i++) //TODO niet in loop berekenen
	{

		//u and v are the vectors within the virtual screen scaled between 0 and 1, so u = px / screenwidth and y = py / screenwidth
		float u = (float)i / (float)screen->width;
		float v = (float)j / (float)screen->height;
		float3 P = view.p1 + u * (view.p2 - view.p1) + v * (view.p3 - view.p1);

		float3 dir = P - view.pos; // vector in the direction you want to shoot your ray
		float3 D = dir / length(dir); //normalize it

		Ray ray = Ray(view.pos, D);

		float3 intersectionColor = Trace(ray);

		screen->pixels[i + j * screen->width] = FloatToIntColor(intersectionColor);
	}
}