#include "BVH.h"

BVH::BVH()
{
	pLeft = nullptr;
	pRight = nullptr;
}

tuple<Vector3f, float> BVH::intersect(Ray & ray, BVH * & bvh)
{
	if(bvh == nullptr)
		return make_tuple(Vector3f(FLT_MAX, FLT_MAX, FLT_MAX), FLT_MAX);
	
	if(bvh->bBox.isIntersected(ray))
	{
		Triangle *pTri = bvh->bBox.pTri;

		if(pTri == nullptr)
		{
			tuple<Vector3f, float> valLeft = intersect(ray, bvh->pLeft);
			tuple<Vector3f, float> valRight = intersect(ray, bvh->pRight);

			if(get<1>(valLeft) < get<1>(valRight))
				return valLeft;
			else
				return valRight;
		}
		else
			return pTri->rayTriangleIntersection(ray);
	}

	return make_tuple(Vector3f(FLT_MAX, FLT_MAX, FLT_MAX), FLT_MAX);
}

void BVH::constructSubtree(vector<Triangle *> & tris, vector<int> & triIdxs, int dim, int begin, int end, BVH **node)
{
	if(begin > end)
		return;

	Vector3f P_min = tris[triIdxs[begin]]->getMinCoords();
	Vector3f P_max = tris[triIdxs[begin]]->getMaxCoords();
	Vector3f P_cog = tris[triIdxs[begin]]->getCenter();

	*node = new BVH();

	for(int i = begin + 1 ; i <= end ; i++)
	{
		Vector3f P_min_ = tris[triIdxs[i]]->getMinCoords();
		Vector3f P_max_ = tris[triIdxs[i]]->getMaxCoords();
		Vector3f P_cog_ = tris[triIdxs[i]]->getCenter();

		P_min(0) = fmin(P_min_(0), P_min(0));
		P_min(1) = fmin(P_min_(1), P_min(1));
		P_min(2) = fmin(P_min_(2), P_min(2));
		P_max(0) = fmax(P_max_(0), P_max(0));
		P_max(1) = fmax(P_max_(1), P_max(1));
		P_max(2) = fmax(P_max_(2), P_max(2));
		P_cog(0) += P_cog_(0);
		P_cog(1) += P_cog_(1);
		P_cog(2) += P_cog_(2);
	}

	(*node)->bBox.P_min = P_min;
	(*node)->bBox.P_max = P_max;

	if(begin == end)
	{
		(*node)->bBox.pTri = tris[triIdxs[begin]];

		return;
	}

	(*node)->bBox.pTri = nullptr;

	float pivot;

	if(dim == 0)
		pivot = P_cog(0) / (end - begin + 1);
	else if(dim == 1)
		pivot = P_cog(1) / (end - begin + 1);
	else
		pivot = P_cog(2) / (end - begin + 1);

	int midIdx;

	if(dim == 0)
		for(int i = begin, j = end ; i <= j ; midIdx = j)
		{
			Vector3f C1 = tris[triIdxs[i]]->getCenter();
			Vector3f C2 = tris[triIdxs[j]]->getCenter();

			if(C1(0) < pivot && C2(0) < pivot)
				i++;
			else if(C1(0) < pivot && C2(0) > pivot)
			{
				i++;
				j--;
			}
			else if(C1(0) > pivot && C2(0) < pivot)
			{
				swap(&triIdxs[i], &triIdxs[j]);
				i++;
				j--;
			}
			else
				j--;
		}
	else if(dim == 1)
		for(int i = begin, j = end ; i <= j ; midIdx = j)
		{
			Vector3f C1 = tris[triIdxs[i]]->getCenter();
			Vector3f C2 = tris[triIdxs[j]]->getCenter();

			if(C1(1) < pivot && C2(1) < pivot)
				i++;
			else if(C1(1) < pivot && C2(1) > pivot)
			{
				i++;
				j--;
			}
			else if(C1(1) > pivot && C2(1) < pivot)
			{
				swap(&triIdxs[i], &triIdxs[j]);
				i++;
				j--;
			}
			else
				j--;
		}
	else
		for(int i = begin, j = end ; i <= j ; midIdx = j)
		{
			Vector3f C1 = tris[triIdxs[i]]->getCenter();
			Vector3f C2 = tris[triIdxs[j]]->getCenter();

			if(C1(2) < pivot && C2(2) < pivot)
				i++;
			else if(C1(2) < pivot && C2(2) > pivot)
			{
				i++;
				j--;
			}
			else if(C1(2) > pivot && C2(2) < pivot)
			{
				swap(&triIdxs[i], &triIdxs[j]);
				i++;
				j--;
			}
			else
				j--;
		}

	if(midIdx == end || midIdx + 1 == begin)
		midIdx = (midIdx + end) / 2;

	dim = (dim + 1) % 3;

	constructSubtree(tris, triIdxs, dim, begin, midIdx, &((*node)->pLeft));
	constructSubtree(tris, triIdxs, dim, midIdx + 1, end, &((*node)->pRight));
}

void BVH::swap(int *idx1, int *idx2)
{
	int *temp = idx1;
	idx1 = idx2;
	idx2 = temp;
}
