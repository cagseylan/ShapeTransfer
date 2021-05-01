#ifndef _MESH_H_
#define _MESH_H_

#include <algorithm>
#include <chrono>
#include <cfloat>
#include <cstdio>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <iostream>
#include <queue>
#include <thread>
#include <tuple>
#include <vector>

#include "Vertex.h"
#include "Triangle.h"
#include "Edge.h"
#include "Skeleton.h"
#include "Ray.h"
#include "BVH.h"
#include "WunderSVD3x3.h"

using namespace Eigen;
using namespace std;

#define SD_ARAP

class Mesh
{
public:
	Mesh();
	Mesh(const char *meshFilePath);	// Mesh construction by reading .off file. 
	Mesh(const char *meshFilePath, const Skeleton *skel);	// Constucts augmented mesh.  
	Mesh(Mesh *sMesh, const Skeleton *sSkel, const Skeleton *tSkel);	// Mesh construction by mesh transfer.

	tuple<Mesh *, Skeleton *> approximate2gtMesh(Mesh *augMesh, const Skeleton *embSkel);
	Mesh *approximate2gtMeshUnconstrained(Mesh *augMesh);
	void saveAsOFF(const char *filePath) const;	// Save the mesh as .off file.
	void savePartitionAsOFF(const char *filePath, int edgeId, int *partitions) const;
	void getPlotPrimitives(MatrixXd & meshVertices, MatrixXi & meshFaces, MatrixXd & faceColors) const;
	void getVerticesFaces(MatrixXd & meshVertices, MatrixXi & meshFaces) const;
	void getDifferenceColorMap(Mesh *approx, MatrixXd & vertexColors);
	unsigned int getNumOfMeshVertices(void) const;
	unsigned int getNumOfMeshEdges(void) const;
	void getSupportVertices(MatrixXd & supportVertices) const;
	vector<tuple<int, int>> getMeshEdges(void) const;
	vector<tuple<int, int, float>> getSupportEdges(void) const;
	vector<float> getEdgeLengths(void) const;
	vector<float> getNonRigidity(void) const;
	vector<float> getEdgeWeights(void) const;
	void setEdgeWeights(vector<float> weights);
	float computeMeshVolume(void);
	vector<vector<float>> computeAPSP(void);
	void postProcess(void);

private:
	int *partitions;
	unsigned long N;	// Number of vertices in the initial mesh. 
	unsigned long M;	// Number of edges in the intial mesh. 
	bool *freeVertices;
	float ratio = 0.5;	// 0.25 for human and horse. 0.48 for wolf. 
	float avgEdgeLength = 0.0;
	float *parametricCoords;
	float meshArea = 0.0;
	vector<int> triIdxs;	// Just for BVH. 
	vector<vector<int>> allSupports;	// i'th skeleton vertex supports allSupports[i][j]'th mesh vertex. 
	vector<float> nonRigidityErrors;	// Nonrigidity per vertex i (Just for transferred meshes). 
	vector<Vertex *> vertices;
	vector<Triangle *> tris;
	vector<Edge *> edges;
	BVH *bvh;

	void copyFromMesh(Mesh *mesh);
	float computeSkelError(Mesh *augMesh, vector<vector<int>> & boneSamples, Skeleton *optSkel);
	vector<vector<int>> assignSamples2Skel(const Skeleton *embSkel);
	vector<int> computeSphereNeighbourhood(Vector3f P, float radius);
	vector<tuple<Vector3f, bool>> computeGtSkelSamples(Mesh *augMesh);
	tuple<Matrix3f, Vector3f> computeRigidTransform(vector<Vector3f> & P, vector<Vector3f> & Q);
	void constructAugmentedMeshV1(const Skeleton *skel);
	void constructAugmentedMeshV2(const Skeleton *skel);
	void constructAugmentedMeshContour(const Skeleton *skel);
	void constructMeshForARAP(const Skeleton *skel);
	void constructAugmentedMeshZhang(const Skeleton *skel);
	void moveHandles(Mesh *sMesh, const Skeleton *sSkel, const Skeleton *tSkel);
	void solveWithARAP(Mesh *sMesh);
	void computeEdgeWeights(void);
	void computeVertexWeights(void);
	Vector3f rayMeshIntersection(Ray & ray) const;
	int nearestVertex(Vector3f P) const;
	void computeEdgeLengths(void);
	void setVerticeNormals(void);
	float computeAxisAlignedAngle(int childEdge, int parentEdge, const Skeleton *sSkel, Vector3f *handles) const;
	void initializeHandles(const Skeleton *sSkel, Vector3f *handles) const;
	float computeOptRAAForTrunk(int trunkIdx, const Skeleton *sSkel, const Skeleton *tSkel) const;
	float computeE_phi(vector<Vector3f> setSrc, vector<int> neighs, const Skeleton *tSkel) const;
	Vector3f rotateAroundEdge(Vector3f P, float theta, int edgeId, const Skeleton *tSkel) const;
	Vector3f rotate(Vector3f P, float theta, Vector3f r) const;
	void parseOFFFile(const char *meshFilePath);
	void correspondMesh(const Skeleton *skel) const;
	tuple<int, float> projectPoint(int meshVertId, const Skeleton *skel) const;
	Vector3f getCartesianCoords(const Skeleton *skel, int skelEdgeId, float t) const;
	Vector3f applyEdgeTranslation(Vector3f P, int edgeId, const Skeleton *sSkel, const Skeleton *tSkel) const;
	Vector3f applyEdgeAlignment(Vector3f P, int edgeId, const Skeleton *sSkel, const Skeleton *tSkel) const;
	Vector3f applyScaling(Vector3f P, int edgeId, float paramCoord, const Skeleton *sSkel, const Skeleton *tSkel) const;
	void printBorderEdgeLengthDiff(const Mesh *sMesh, int *partitions);
	void markFreeVertices(float ratio, bool *freeVertices, const Skeleton *skel) const;
	float computeProjEdgeWeight(int meshVertId) const;
	void saveNonRigidPartsAsOFF(const char *filePath, bool *freeVertices) const;
	float computeCotan(int edgeId) const;
	void computeOptimumRotation(const Mesh *sMesh, int idx, SparseMatrix<float> & W, vector<vector<int>> & neighbourhoodInfo, Matrix3f & R) const;
	void computeOptimumRotation_(const Mesh *sMesh, SparseMatrix<float> & W, vector<vector<int>> & neighbourhoodInfo, vector<Matrix3f> & R, int b, int e) const;
	Matrix3f computeOptimumSmoothRotation(const Mesh *sMesh, int idx, SparseMatrix<float> & W, vector<Matrix3f> & Rs) const;
	void updatePositions(const Mesh *sMesh, vector<Matrix3f> & R, SparseLU<SparseMatrix<float>, COLAMDOrdering<int>> & solver, SparseMatrix<float> & W, SparseMatrix<float> & L, bool *freeVertices) const;
	void solveForPositions(VectorXf & b, SparseLU<SparseMatrix<float>, COLAMDOrdering<int>> & solver, VectorXf & x) const;
	float computeError(const Mesh *sMesh, vector<Matrix3f> & R, SparseMatrix<float> & W, bool onlyMesh) const;
	float computePerCellError(const Mesh *sMesh, int vert, Matrix3f & R, SparseMatrix<float> & W, bool onlyMesh) const;
	void taubinMeshFairingOneStep(void);
	void computeInitialSolution(Mesh *sMesh, const Skeleton *sSkel, const Skeleton *tSkel);
	Matrix4f computeTransformationMatrix(const Skeleton *sSkel, const Skeleton *tSkel, float phi, int edgeId);
	void loadCoeffs(const char *filePath, vector<vector<float>> & coeffs);
};

#endif
