#include <chrono>
#include <cstdlib>
#include <Eigen/Dense>
#include <igl/opengl/glfw/Viewer.h>
#include <tuple>
#include <vector>

#include "Mesh.h"
#include "Skeleton.h"

using namespace Eigen;
using namespace std;

void saveSupportEdges(MatrixXd & meshVertices, MatrixXd & supportVertices, vector<tuple<int, int, float>> & supportEdges, const char *fileName)
{
	FILE *pFile = fopen(fileName, "w");

	for(int i = 0 ; i < supportEdges.size() ; i++)

		if(get<2>(supportEdges[i]) != 0.0)
		{
			int p1idx = get<0>(supportEdges[i]) - meshVertices.rows();
			int p2idx = get<1>(supportEdges[i]);

			float P1x = supportVertices.row(p1idx)(0);
			float P1y = supportVertices.row(p1idx)(1);
			float P1z = supportVertices.row(p1idx)(2);
			float P2x = meshVertices.row(p2idx)(0);
			float P2y = meshVertices.row(p2idx)(1);
			float P2z = meshVertices.row(p2idx)(2);

			fprintf(pFile, "%f %f %f %f %f %f\n", P1x, P1y, P1z, P2x, P2y, P2z);
		}

	fclose(pFile);
}

int main(int argc, char **argv)
{
	const Skeleton *sourceSkel = new Skeleton(argv[1]);	// Embedded skeleton. 
	printf("Source skeleton loaded.\n");

	const Skeleton *targetSkel = new Skeleton(argv[2]);	// Target skeleton (into which the source mesh will be transferred)
	printf("Target skeleton loaded.\n");

	printf("Comparison of source and target skeleton edge lengths:\n");
	for(int i = 0 ; i < sourceSkel->edges.size() ; i++)
		printf("%d - %f %f\n", i, sourceSkel->edges[i]->length, targetSkel->edges[i]->length);

	//Mesh *gtMesh = new Mesh(argv[4]);
	Mesh *sourceMesh = new Mesh(argv[3], sourceSkel);	// Source mesh, argv[1] is the emb. skeleton of the source mesh.
	Mesh *targetMesh;

	bool showOnlySource = false;
	bool showGroundTruth = false;
	unsigned long N;
	unsigned long M;
	MatrixXd meshVertices;
	MatrixXd skeletonVertices;
	MatrixXd supportVertices;
	MatrixXd faceColors;
	MatrixXd vertexColors;
	MatrixXi meshFaces;
	vector<tuple<int, int>> skeletonEdges;
	vector<tuple<int, int, float>> supportEdges;

	if(showGroundTruth)
	{
		//gtMesh->getVerticesFaces(meshVertices, meshFaces);
	}
	else if(showOnlySource)
	{
		sourceMesh->getPlotPrimitives(meshVertices, meshFaces, faceColors);

		sourceSkel->getSkeletonVertices(skeletonVertices);
		skeletonEdges = sourceSkel->getSkeletonEdges();

		N = sourceMesh->getNumOfMeshVertices();
		M = sourceMesh->getNumOfMeshEdges();

		sourceMesh->getSupportVertices(supportVertices);
		supportEdges = sourceMesh->getSupportEdges();

		// saveSupportEdges(meshVertices, supportVertices, supportEdges, argv[4]);
	}
	else
	{
		targetMesh = new Mesh(sourceMesh, sourceSkel, targetSkel);

		//targetMesh->postProcess();
		
		//targetMesh = gtMesh->approximate2gtMeshUnconstrained(sourceMesh);
		
		//tuple<Mesh *, Skeleton *> approxRecord = gtMesh->approximate2gtMesh(sourceMesh, sourceSkel);
		//targetMesh = get<0>(approxRecord);
		//Skeleton *optSkel = get<1>(approxRecord);
		//gtMesh->getDifferenceColorMap(targetMesh, vertexColors);

		//targetMesh->saveAsOFF("../tr_reg_008_approxMesh.off");
		//optSkel->saveAsCG("../tr_reg_008_optSkel.cg");

		//targetMesh->saveAsOFF("../tr_reg_003_reconstructed.off");
		
		// --------------------------------------- Volume computations --------------------------------------------
		float initialVol = sourceMesh->computeMeshVolume();
		float resultVol = targetMesh->computeMeshVolume();
		float change = 100.0 * (resultVol - initialVol) / initialVol;

		printf("Initial Volume: %f Result Volume: %f Change: %f\n", initialVol, resultVol, change);
		// --------------------------------------------------------------------------------------------------------

		// ---------------------------------- Isometric error computations ----------------------------------------
		/*vector<vector<float>> initialAPSP = sourceMesh->computeAPSP();
		vector<vector<float>> resultAPSP = targetMesh->computeAPSP();
		float isometricErr = 0.0;

		for(int i = 0 ; i < initialAPSP.size() ; i++)
			for(int j = 0 ; j < initialAPSP[i].size() ; j++)
			{
				float diff = resultAPSP[i][j] - initialAPSP[i][j];

				isometricErr += diff * diff;
			}

		isometricErr = sqrt(isometricErr);

		printf("Isometric error: %f\n", isometricErr);*/
		// --------------------------------------------------------------------------------------------------------

		// ---------------------------- Primitives required to plot the resulting mesh ----------------------------
		targetMesh->getPlotPrimitives(meshVertices, meshFaces, faceColors);
		// --------------------------------------------------------------------------------------------------------

		// --------------------------- Primitives required to plot the curve skeleton -----------------------------
		targetSkel->getSkeletonVertices(skeletonVertices);
		skeletonEdges = targetSkel->getSkeletonEdges();
		// --------------------------------------------------------------------------------------------------------
	
		// ---------------------------- Primitives required to plot the support edges -----------------------------
		//
		N = targetMesh->getNumOfMeshVertices();
		M = targetMesh->getNumOfMeshEdges();

		targetMesh->getSupportVertices(supportVertices);
		supportEdges = targetMesh->getSupportEdges();
		// --------------------------------------------------------------------------------------------------------

		// -------------------- Computation of the distortion per 3d patch after transfer -------------------------
		/*MatrixXd nonRigidityColors(N, 3);
		vector<float> nonRigidityErrors = targetMesh->getNonRigidity();
	
		vector<float>::iterator maxNonRigidityIter = max_element(nonRigidityErrors.begin(), nonRigidityErrors.end());
		float maxNonRigidity = *maxNonRigidityIter;

		for(int i = 0 ; i < N ; i++)
		{
			if(maxNonRigidity == 0.0)
				maxNonRigidity = 0.000001;

			nonRigidityColors(i, 0) = nonRigidityErrors[i] / maxNonRigidity;
			nonRigidityColors(i, 1) = 0;
			nonRigidityColors(i, 2) = 1 - nonRigidityErrors[i] / maxNonRigidity;
		}*/
		// --------------------------------------------------------------------------------------------------------

		// ---------------------------------- Computation of the difference in edge lengths -----------------------
		/*vector<tuple<int, int>> meshEdges = sourceMesh->getMeshEdges();

		vector<float> edgeLengthDiffs;
		vector<float> srcEdgeLengths = sourceMesh->getEdgeLengths();
		vector<float> trgEdgeLengths = targetMesh->getEdgeLengths();

		for(int i = 0 ; i < srcEdgeLengths.size() ; i++)
			edgeLengthDiffs.push_back(fabs(srcEdgeLengths[i] - trgEdgeLengths[i]));

		vector<float>::iterator maxEdgeDiffIter = max_element(edgeLengthDiffs.begin(), edgeLengthDiffs.end() - M);
		float maxEdgeDiff = *maxEdgeDiffIter;*/
		// --------------------------------------------------------------------------------------------------------
	}

	igl::opengl::glfw::Viewer viewer;

	viewer.core().background_color.setOnes();

	// Plot the transferred mesh. 
	viewer.data().set_mesh(meshVertices, meshFaces);
	viewer.data().set_colors(faceColors);

	// Uncomment below to plot the non-rigidity error per vertex. 
	//viewer.data().set_colors(nonRigidityColors);

	// Uncomment below to plot the edge length differences. 
	/*for(int i = 0 ; i < meshEdges.size() ; i++)
	{
		int u = get<0>(meshEdges[i]);
		int v = get<1>(meshEdges[i]);

		float redColor = (edgeLengthDiffs[2 * i]) / maxEdgeDiff;

		viewer.data().add_edges(meshVertices.row(u), meshVertices.row(v), RowVector3d(redColor, 0, 0));
	}*/

	//Uncomment below to plot the support edges. 
	int numOfSupports = 0;

	for(int i = 0 ; i < supportEdges.size() ; i++)
	{
		float w = get<2>(supportEdges[i]);

		if(w != 0.0)
		{
			numOfSupports++;

			int u = get<0>(supportEdges[i]) - N;
			int v = get<1>(supportEdges[i]);

			viewer.data().add_edges(supportVertices.row(u), meshVertices.row(v), RowVector3d(0, 1, 0));
		}
	}

	printf("Num. of mesh vertices: %lu\n", N);
	printf("Num. of mesh edges: %lu\n", M / 2);
	printf("Num. of supports: %d\n", numOfSupports / 2);

	// Uncomment below to plot the skeleton edges. 
	for(int i = 0 ; i < skeletonEdges.size() ; i++)
	{
		int u = get<0>(skeletonEdges[i]);
		int v = get<1>(skeletonEdges[i]);

		viewer.data().add_edges(skeletonVertices.row(u), skeletonVertices.row(v), RowVector3d(1, 0, 0));
	}

	// Gouraud shading
	// viewer.data().set_face_based(true); 

	viewer.launch();

	return 0;
}
