#include "Triangle.h"

Triangle::Triangle(int id, int v1, int v2, int v3)
	: id(id), v1(v1), v2(v2), v3(v3)
{}

/*
	Given id's of two vertices, it returns the id of the third vertex. 
*/
int Triangle::getThirdVertex(int vert1, int vert2)
{
	return v1 + v2 + v3 - vert1 - vert2;
}

/*
	Sets normal of the triangle as a unit vector. 
	Assumed that v1, v2, v3 are in CCW order. 
*/
void Triangle::setNormal(Vector3f V1, Vector3f V2, Vector3f V3)
{
	normal = (V2 - V1).cross(V3 - V1);
	normal.normalize();
}

tuple<Vector3f, float> Triangle::rayTriangleIntersection(Ray & ray)
{
	Vector3f o = ray.o;
	Vector3f d = ray.d;
	Vector3f P1 = (*vertices)[this->v1]->coords;
	Vector3f P2 = (*vertices)[this->v2]->coords;
	Vector3f P3 = (*vertices)[this->v3]->coords;
	float xa_m_xe = P1(0) - o(0);	
	float ya_m_ye = P1(1) - o(1);
	float za_m_ze = P1(2) - o(2);
	float det_A = this->xa_m_xb * (this->ya_m_yc * d(2) - d(1) * this->za_m_zc) - this->xa_m_xc * (this->ya_m_yb * d(2) - d(1) * this->za_m_zb) + d(0) * this->det_A_term;
	const float intTestEps = 1e-4;

	float t = (this->xa_m_xb * (this->ya_m_yc * za_m_ze - ya_m_ye * this->za_m_zc) - this->xa_m_xc * (this->ya_m_yb * za_m_ze - ya_m_ye * this->za_m_zb) + xa_m_xe * this->det_A_term) / det_A;

	if(t < -intTestEps)
		return make_tuple(Vector3f(FLT_MAX, FLT_MAX, FLT_MAX), FLT_MAX);

	float gamma = (this->xa_m_xb * (ya_m_ye * d(2) - d(1) * za_m_ze) - xa_m_xe * (this->ya_m_yb * d(2) - d(1) * this->za_m_zb) + d(0) * (this->ya_m_yb * za_m_ze - ya_m_ye * this->za_m_zb)) / det_A;

	if(gamma < -intTestEps || gamma > 1 + intTestEps)
		return make_tuple(Vector3f(FLT_MAX, FLT_MAX, FLT_MAX), FLT_MAX);

	float beta = (xa_m_xe * (this->ya_m_yc * d(2) - d(1) * this->za_m_zc) - this->xa_m_xc * (ya_m_ye * d(2) - d(1) * za_m_ze) + d(0) * (ya_m_ye * this->za_m_zc - this->ya_m_yc * za_m_ze)) / det_A;

	if(beta < -intTestEps || beta > 1 - gamma + intTestEps)
		return make_tuple(Vector3f(FLT_MAX, FLT_MAX, FLT_MAX), FLT_MAX);

	return make_tuple(o + d * t, t);
}

Vector3f Triangle::getMinCoords(void)
{
	Vector3f P_min;

	Vector3f a = (*vertices)[this->v1]->coords;
	Vector3f b = (*vertices)[this->v2]->coords;
	Vector3f c = (*vertices)[this->v3]->coords;

	P_min(0) = fmin(fmin(a(0), b(0)), c(0));
	P_min(1) = fmin(fmin(a(1), b(1)), c(1));
	P_min(2) = fmin(fmin(a(2), b(2)), c(2));

	return P_min;
}

Vector3f Triangle::getMaxCoords(void)
{
	Vector3f P_max;

	Vector3f a = (*vertices)[this->v1]->coords;
	Vector3f b = (*vertices)[this->v2]->coords;
	Vector3f c = (*vertices)[this->v3]->coords;

	P_max(0) = fmax(fmax(a(0), b(0)), c(0));
	P_max(1) = fmax(fmax(a(1), b(1)), c(1));
	P_max(2) = fmax(fmax(a(2), b(2)), c(2));

	return P_max;
}

Vector3f Triangle::getCenter(void)
{
	Vector3f a = (*vertices)[this->v1]->coords;
	Vector3f b = (*vertices)[this->v2]->coords;
	Vector3f c = (*vertices)[this->v3]->coords;

	return (a + b + c) / 3.0;
}
