#ifndef _BOX_H_
#define _BOX_H_

#include <Eigen/Dense>

#include "Ray.h"
#include "Triangle.h"

using namespace Eigen;

class Box
{
public:
	bool isIntersected(Ray & ray);

	Vector3f P_min;
	Vector3f P_max;
	Triangle *pTri = nullptr;
};

#endif
