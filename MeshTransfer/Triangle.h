#ifndef _TRIANGLE_H_
#define _TRIANGLE_H_

#include <cfloat>
#include <Eigen/Dense>
#include <vector>

#include "Ray.h"
#include "Vertex.h"

using namespace Eigen;
using namespace std;

class Triangle
{
public:
	int id;
	int v1;
	int v2;
	int v3;
	Vector3f normal;
	vector<Vertex *> *vertices;

	// Below are only for computing min. dist. to a point and ray-triangle intersection
	Vector3f U1;
	Vector3f U2;
	Vector3f U3;
	float xa_m_xb;
	float xa_m_xc;
	float ya_m_yb;
	float ya_m_yc;
	float za_m_zb;
	float za_m_zc;
	float det_A_term;

	Triangle(int id, int v1, int v2, int v3);
	int getThirdVertex(int vert1, int vert2);
	void setNormal(Vector3f V1, Vector3f V2, Vector3f v3);
	tuple<Vector3f, float> rayTriangleIntersection(Ray & ray);
	Vector3f getMinCoords(void);
	Vector3f getMaxCoords(void);
	Vector3f getCenter(void);
};

#endif
