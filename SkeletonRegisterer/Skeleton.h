#ifndef _SKELETON_H_
#define _SKELETON_H_

#include <cfloat>
#include <cstdio>
#include <Eigen/Dense>
#include <queue>
#include <tuple>
#include <vector>

#include "Edge.h"
#include "Vertex.h"

using namespace Eigen;

class Skeleton
{
public:
	int shortestEdgeIdx = -1;	// Id of the shortest edge in the skeleton. 
	double totalLength = 0.0;	// Total length of the skeleton. 
	vector<vector<int>> branchInfo;
	vector<vector<int>> segments;
	vector<Edge *> edges;
	vector<Vertex *> vertices;

	Skeleton(void);
	Skeleton(const char *fileName);
	Skeleton(const Skeleton *embSkel, const Skeleton *curveSkel, vector<int> pairs);
	void saveAsOFF(const char *fileName);	// Saves the vertices in OFF format. TODO: Subject to be removed.
	void saveAsCG(const char *fileName);	// Saves the skeleton in CG format into the path fileName.
	double computeDissimilarity(const Skeleton *skel) const;	// Computes similarity between two curve-skeletons.
	double computeHausdorffDistance(const Skeleton *skel) const;	// Computes Hausdorff distance between two curve-skeletons. 
	vector<int> computeShortestPath_BFS(int source, int dest, double *distanceInfo) const;	// Computes shortest path between source and dest.
	void getSkeletonVertices(MatrixXd & skeletonVertices) const;
	vector<tuple<int, int>> getSkeletonEdges(void) const;

private:
	void computeBranchInfo(void);	// Fills branchInfo vector.
	void computeSegments(void);	// Computes skeleton segments.  
	double computeDistance2Skeleton(const Vector3f & point) const;	// Computes distance between the skeleton and a point.
	double computeDistance2LineSeg(const Vector3f & U, const Vector3f & V, const Vector3f & P) const; // Computes distance btw. a line segment and a point. 
	double computeAvgDistance(double stepSize, const Skeleton *skel) const;	//  Computes average distance between skel and this skeleton.
	tuple<double, int> computeDistInfoOfSegment(double stepSize, vector<int> segment, const Skeleton * & skel) const; 
	double supInt(double stepSize, const Skeleton *skel) const;	// Computes Sup Int terms in Hausdorff distance with stepSize.
	vector<Vector3f> sampleSegment(double stepSize, vector<int> segment) const; // Samples a skeleton segment (which is an arc) with a spacing of stepSize.
};

#endif
