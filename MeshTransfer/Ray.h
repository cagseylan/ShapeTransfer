#ifndef _RAY_H_
#define _RAY_H_

#include <Eigen/Dense>

using namespace Eigen;

class Ray
{
public:
	Vector3f o;	// Origin of the ray
	Vector3f d;	// Direction of the ray

	Ray(Vector3f o, Vector3f d);
};

#endif
