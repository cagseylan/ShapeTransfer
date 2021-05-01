#ifndef _SKELETON_H_
#define _SKELETON_H_

#include <cstdio>
#include <Eigen/Dense>
#include <vector>

#include "Edge.h"
#include "Vertex.h"

using namespace Eigen;
using namespace std;

class Skeleton
{
public:
	float avgEdgeLength = 0.0;
	vector<Vertex *> vertices;
	vector<Edge *> edges;

	Skeleton();
	Skeleton(const char *filePath);	// Skeleton construction by reading .cg file. 
	void saveAsCG(const char *fileName);	// Saves the skeleton in CG format into the path fileName.
	int findTrunk(void) const;
	void getSkeletonVertices(MatrixXd & skeletonVertices) const;
	vector<tuple<int, int>> getSkeletonEdges(void) const;
	vector<tuple<Vector3f, Vector3f>> sampleSkeleton(float stepSize) const;
	vector<tuple<Vector3f, vector<Vector3f>>> computeRays(float stepSize, int raysPerSample) const;
	double computeDistance2LineSeg(const Vector3f & U, const Vector3f & V, const Vector3f & P) const; // Computes distance btw. a line segment and a point.

private:
	float edgeLength(int edgeId); // Returns length of the edge with edgeId.
	vector<vector<int>> computeSegments(void) const;	// Computes skeleton segments. 
	vector<tuple<Vector3f, Vector3f>> sampleSegment(vector<int> segment, float stepSize) const;
	vector<Vector3f> computeDirections(Vector3f origin, Vector3f edgeDirection, int raysPerSample) const;
};

#endif
