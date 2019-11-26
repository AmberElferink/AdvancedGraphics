#include "core_settings.h"

// adapted from Möller–Trumbore intersection algorithm: https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
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
		intersection.material = scene.matList[triangle.material];
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
				if (length(intersection.point - ray.O) <  length(ray.O - light.position)) //Between the light and the origin, not after
					return true;
		}
		return false; //false if no intersections are found
	}
}

/*Method that shoots a shadow ray and checks whether there are objects between the current intersection point and a light source*/
bool Raytracer::viewLight( Intersection intersection, const Light &light, float3 &lightVector )
{
	float3 dir = intersection.point - light.position; //vector between light and intersection point
	float3 D = dir / length( dir ); //normalized vector

	lightVector = D;

	Ray shadowRay = Ray( intersection.point + intersection.norm * 0.001f, D ); //shadow ray from origin to light point

	Intersection closest;

	if ( IsOccluded( shadowRay, light ) ) 
		return false; //cannot see light source
	else
		return true; //no objects that obstruct view of light source
}

uint Raytracer::FloatToIntColor( float3 floatColor )
{
	return ( (int)( floatColor.x * 255.0f ) << 16 + (int)( floatColor.y * 255.0f ) << 8 + (int)( floatColor.z * 255.0f ) );
}

//-----------------------------------------------------
// shoot a ray through point p to intersect the scene.
//  p3 |------------------| 
//     |        .         |
//     |        p         | screenheight
//  p1 |------------------| p2
//   screenwidth
//-----------------------------------------------------
void Raytracer::rayTrace(Bitmap *screen, const ViewPyramid &view)
{
	Timer t;
	for ( int j = 0; j < screen->height; j++ )
	{
		for (int i = 0; i < screen->width; i++) //TODO niet in loop berekenen
		{

			//u and v are the vectors within the virtual screen scaled between 0 and 1, so u = px / screenwidth and y = py / screenwidth
			float u = (float)i / (float)screen->width; 
			float v = (float)j / (float)screen->height;
			float3 P = view.p1 + u * ( view.p2 - view.p1 ) + v * ( view.p3 - view.p1 );

			float3 dir = P - view.pos; // vector in the direction you want to shoot your ray
			float3 D = dir / length( dir ); //normalize it

			Ray ray = Ray( view.pos, D ); 

			Intersection closest; //this will be your closest intersection of which you want to know the color
			closest.t = 10e30;
			closest.material = new Material( make_float3( 0, 0, 0 ) ); //default black (background)

			//Find closest intersection point for all meshes
			for ( Mesh &mesh : scene.meshList )
			{
				t.reset();
				int triangleCount = mesh.vcount / 3;
				for ( int i = 0; i < triangleCount; i++ ) //find the closest triangle intersection for all triangles
				{
					Intersection intersection;
					if (Intersect( ray, mesh.triangles[i], intersection ) )
					{
						if ( intersection.t < closest.t ) //update closest intersection point only if distance to intersection point is shorter
							closest = intersection; 
					}
				}
			}

			float3 intersectionColor = make_float3(0,0,0);

			/*Check for all lights whether they can be seen from the current intersection point*/
			for (Light &light : scene.lightList)
			{
				float3 lightVector;
				if ( viewLight( closest, light, lightVector ) )
					intersectionColor += closest.material->diffuse * light.radiance * dot( closest.norm, lightVector ); //If light source can be seen, multiply color with current pixel color
			}

			screen->pixels[i + j * screen->width] = FloatToIntColor( intersectionColor );
		}
	}
}