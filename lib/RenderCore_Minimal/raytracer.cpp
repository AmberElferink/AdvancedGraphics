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