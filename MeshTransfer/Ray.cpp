#include "Ray.h"

Ray::Ray(Vector3f o, Vector3f d)
{
	this->o = o;
	this->d = d.normalized();
}
