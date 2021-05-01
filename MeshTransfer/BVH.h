#ifndef _BVH_H_
#define _BVH_H_

#include <cfloat>
#include <Eigen/Dense>
#include <tuple>

#include "Triangle.h"
#include "Box.h"

class BVH
{
public:
	BVH();
	tuple<Vector3f, float> intersect(Ray & ray, BVH * & bvh);
	void constructSubtree(vector<Triangle *> & tris, vector<int> & triIdxs, int dim, int begin, int end, BVH **node);
	void swap(int *idx1, int *idx);

private:
	BVH *pLeft;
	BVH *pRight;
	Box bBox;
};

#endif
