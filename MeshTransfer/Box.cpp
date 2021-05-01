#include "Box.h"

bool Box::isIntersected(Ray & ray)
{
	float tmin;
	float tmax;
	float tymin;
	float tymax; 
	float tzmin;
	float tzmax; 

	Vector3f dir = ray.d;
	Vector3f org = ray.o;

	if(dir(0) >= 0)
	{
		tmin = (P_min(0) - org(0)) / dir(0);
		tmax = (P_max(0) - org(0)) / dir(0);
	}
	else
	{
		tmin = (P_max(0) - org(0)) / dir(0);
		tmax = (P_min(0) - org(0)) / dir(0);
	}

	if(dir(1) >= 0)
	{
		tymin = (P_min(1) - org(1)) / dir(1);
		tymax = (P_max(1) - org(1)) / dir(1);
	}
	else
	{
		tymin = (P_max(1) - org(1)) / dir(1);
		tymax = (P_min(1) - org(1)) / dir(1);
	}

	if((tmin > tymax) || (tymin > tmax))
		return false;

	if(tymin > tmin)
		tmin = tymin;

	if(tymax < tmax)
		tmax = tymax;

	if(dir(2) >= 0)
	{
		tzmin = (P_min(2) - org(2)) / dir(2);
		tzmax = (P_max(2) - org(2)) / dir(2);
	}

	else
	{
		tzmin = (P_max(2) - org(2)) / dir(2);
		tzmax = (P_min(2) - org(2)) / dir(2);
	}

	if((tmin > tzmax) || (tzmin > tmax))
		return false;

	if(tzmin > tmin)
		tmin = tzmin;

	if(tzmax < tmax)
		tmax = tzmax;

	return true;
}
