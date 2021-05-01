#include "Mesh.h"

Mesh::Mesh()
{}

/*
	Mesh construction by reading .off file. 
*/
Mesh::Mesh(const char *meshFilePath)
{
	parseOFFFile(meshFilePath);

	computeEdgeLengths();
	computeEdgeWeights();
	computeVertexWeights();

	partitions = (int *) malloc(sizeof(int) * vertices.size());
	parametricCoords = (float *) malloc(sizeof(float) * vertices.size());
}

/*
	Constucts augmented mesh. 
*/
Mesh::Mesh(const char *meshFilePath, const Skeleton *skel)
{
	parseOFFFile(meshFilePath);

	auto start = std::chrono::steady_clock::now();

	printf("Bounding volume hierarchy is being constructed...\n");

	bvh->constructSubtree(tris, triIdxs, 0, 0, tris.size() - 1, &bvh);

	auto end = std::chrono::steady_clock::now();

	std::chrono::duration<double> elapsedSeconds = end - start;

	printf("Bounding volume hierarchy constructed. (%lf seconds)\n", elapsedSeconds.count());

	constructAugmentedMeshZhang(skel);

/*
#ifdef SD_ARAP
	constructAugmentedMeshContour(skel);
#else
	constructMeshForARAP(skel);
#endif
*/

	computeVertexWeights();
}

/* 
	Mesh construction by mesh transfer. 
	sMesh: source mesh. 
	tSkel: target skeleton. 
*/
Mesh::Mesh(Mesh *sMesh, const Skeleton *sSkel, const Skeleton *tSkel)
{
	// ----------------------------- The connection of the transferred mesh is same as the source mesh ----------------------------------
	this->N = sMesh->N;
	this->M = sMesh->M;
	// Number of vertices are same with the source mesh. Deep copy. 
	for(int i = 0 ; i < sMesh->vertices.size() ; i++)
	{
		vertices.push_back(new Vertex(i, sMesh->vertices[i]->coords));

		// Outgoing edge lists are same as the source mesh. Deep copy. 
		for(int j = 0 ; j < sMesh->vertices[i]->outgoingEdges.size() ; j++) 
			vertices[i]->outgoingEdges.push_back(sMesh->vertices[i]->outgoingEdges[j]);
	}

	// Connection of the triangles is same with the source mesh. Deep copy. 
	for(int i = 0 ; i < sMesh->tris.size() ; i++)
		tris.push_back(new Triangle(sMesh->tris[i]->id, sMesh->tris[i]->v1, sMesh->tris[i]->v2, sMesh->tris[i]->v3));

	// Edges are also same with the source mesh. Deep copy.
	for(int i = 0 ; i < sMesh->edges.size() ; i++)
		edges.push_back(new Edge(sMesh->edges[i]->id, sMesh->edges[i]->u, sMesh->edges[i]->v, sMesh->edges[i]->opp, sMesh->edges[i]->faceId));

	// Arrange partitions and parametric coordinates. They will be same with the source mesh. 
	freeVertices = (bool *) malloc(sizeof(bool) * vertices.size());
	partitions = (int *) malloc(sizeof(int) * vertices.size());
	parametricCoords = (float *) malloc(sizeof(float) * vertices.size());

	for(int i = 0 ; i < vertices.size() ; i++)
	{
		freeVertices[i] = sMesh->freeVertices[i];
		partitions[i] = sMesh->partitions[i];
		parametricCoords[i] = sMesh->parametricCoords[i];
	}

	for(int i = 0 ; i < edges.size() ; i++)
		edges[i]->w = sMesh->edges[i]->w;

	ratio = sMesh->ratio;
	// ----------------------------------------------------------------------------------------------------------------------------------

	printf("Moving handles...\n");

	auto start = std::chrono::steady_clock::now();

	moveHandles(sMesh, sSkel, tSkel);	// Move the mesh onto the target skeleton

	auto end = std::chrono::steady_clock::now();

	std::chrono::duration<double> elapsedSeconds = end - start;

	printf("Handles moved. (%lf seconds)\n", elapsedSeconds.count());

	// ----------------------------------------------------------------------------------------------------------------------------------

	//computeInitialSolution(sMesh, sSkel, tSkel);

	solveWithARAP(sMesh);

	computeEdgeLengths();
}

/*
	Method to compute approximate mesh to the ground truth mesh. Computes both the best approximation, with the best target skeleton 
	giving this best approximation. 
	augMesh: The augmented source mesh
	embSkel: Embedded skeleton of the source mesh
*/
tuple<Mesh *, Skeleton *> Mesh::approximate2gtMesh(Mesh *augMesh, const Skeleton *embSkel)
{	
	Mesh *approximateMesh = new Mesh();
	Skeleton *optSkel = new Skeleton();

	vector<tuple<Vector3f, bool>> newSampleRecords = this->computeGtSkelSamples(augMesh);
	vector<vector<int>> boneSamples = augMesh->assignSamples2Skel(embSkel);
	vector<tuple<Matrix3f, Vector3f>> optBoneTransformations;
	vector<int> numMeaningfulSamples(boneSamples.size());

	for(int i = 0 ; i < numMeaningfulSamples.size() ; i++)
		numMeaningfulSamples[i] = 0;

	// Compute optimum bone transformations. 
	for(int i = 0 ; i < boneSamples.size() ; i += 2)
	{
		vector<Vector3f> P;
		vector<Vector3f> Q;

		for(int j = 0 ; j < boneSamples[i].size() ; j++)
		{
			bool isMeaningful = get<1>(newSampleRecords[boneSamples[i][j]]);

			if(isMeaningful)
			{
				P.push_back(augMesh->vertices[boneSamples[i][j] + augMesh->N]->coords);
				Q.push_back(get<0>(newSampleRecords[boneSamples[i][j]]));

				numMeaningfulSamples[i]++;
				numMeaningfulSamples[i + 1]++;
			}
		}

		tuple<Matrix3f, Vector3f> transformation = this->computeRigidTransform(P, Q);

		Matrix3f R = get<0>(transformation);
		Vector3f t = get<1>(transformation);

		optBoneTransformations.push_back(transformation);
		optBoneTransformations.push_back(transformation);
	}

	// optSkel is same with embSkel topologically. 
	for(int i = 0 ; i < embSkel->vertices.size() ; i++)
	{
		optSkel->vertices.push_back(new Vertex(i, embSkel->vertices[i]->coords));
		optSkel->vertices[i]->outgoingEdges = embSkel->vertices[i]->outgoingEdges;
	}

	for(int i = 0 ; i < embSkel->edges.size() ; i++)
		optSkel->edges.push_back(new Edge(i, embSkel->edges[i]->u, embSkel->edges[i]->v, embSkel->edges[i]->opp, -1));

	//  Move the vertices.
	for(int i = 0 ; i < embSkel->vertices.size() ; i++)
	{
		Vector3f vertPos = embSkel->vertices[i]->coords;
		vector<Vector3f> transformedPosList;
		int denom = 0;

		for(int j = 0 ; j < embSkel->vertices[i]->outgoingEdges.size() ; j++)
		{
			int edgeIdx = embSkel->vertices[i]->outgoingEdges[j];

			Matrix3f R = get<0>(optBoneTransformations[edgeIdx]);
			Vector3f t = get<1>(optBoneTransformations[edgeIdx]);

			transformedPosList.push_back(numMeaningfulSamples[edgeIdx] * (R * vertPos + t));

			denom += numMeaningfulSamples[edgeIdx];
		}

		Vector3f optPose(0.0, 0.0, 0.0);

		for(int j = 0 ; j < transformedPosList.size() ; j++)
			optPose += transformedPosList[j];

		optPose /= denom;

		optSkel->vertices[i]->coords = optPose;
	}

	float skelError = this->computeSkelError(augMesh, boneSamples, optSkel);
	printf("Total Skeleton Error: %f\n", skelError);

	approximateMesh->copyFromMesh(augMesh);

	approximateMesh->moveHandles(augMesh, embSkel, optSkel);

	for(int i = 0 ; i < augMesh->N ; i++)
		if(augMesh->freeVertices[i])
			approximateMesh->vertices[i]->coords = augMesh->vertices[i]->coords;
		else
			approximateMesh->vertices[i]->coords = this->vertices[i]->coords;

	// ------------------------------------------------------------------------------
	// Do not stitch the joints. 
	/*approximateMesh->copyFromMesh(augMesh);

	for(int i = 0 ; i < boneSamples.size() ; i += 2)
	{	
		for(int j = 0 ; j < boneSamples[i].size() ; j++)
		{
			bool isMeaningful = get<1>(newSampleRecords[boneSamples[i][j]]);

			if(isMeaningful)
			{
				int vertId = boneSamples[i][j] + augMesh->N;
				Matrix3f R = get<0>(optBoneTransformations[i]);
				Vector3f t = get<1>(optBoneTransformations[i]);

				approximateMesh->vertices[vertId]->coords = R * augMesh->vertices[vertId]->coords + t;
			}
		}
	}

	for(int i = 0 ; i < augMesh->N ; i++)
		if(augMesh->freeVertices[i])
			approximateMesh->vertices[i]->coords = augMesh->vertices[i]->coords;
		else
			approximateMesh->vertices[i]->coords = this->vertices[i]->coords;*/
	// ------------------------------------------------------------------------------

	approximateMesh->solveWithARAP(augMesh);

	approximateMesh->computeEdgeLengths();

	// ------------------------------------------------------------------------------
	for(int k = 0 ; k < 10 ; k++)
	{
		vector<vector<int>> closestVertices(embSkel->vertices.size(), vector<int>());

		for(int i = 0 ; i < augMesh->N ; i++)
		{
			int closestIdx;
			float minDist = FLT_MAX;
			Vector3f meshVertCoors = augMesh->vertices[i]->coords;

			for(int j = 0 ; j < embSkel->vertices.size() ; j++)
			{
				Vector3f skelVertCoords = embSkel->vertices[j]->coords;
				float currDist = (meshVertCoors - skelVertCoords).norm();

				if(currDist < minDist)
				{
					minDist = currDist;
					closestIdx = j;
				}
			}

			closestVertices[closestIdx].push_back(i);
		}

		vector<Vector3f> displacements;

		for(int i = 0 ; i < approximateMesh->N ; i++)
			displacements.push_back(this->vertices[i]->coords - approximateMesh->vertices[i]->coords);

		for(int i = 0 ; i < optSkel->vertices.size() ; i++)
		{
			Vector3f avgDispVect(0.0, 0.0, 0.0);

			for(int j = 0 ; j < closestVertices[i].size() ; j++)
				avgDispVect += displacements[closestVertices[i][j]];

			avgDispVect /= closestVertices[i].size();

			optSkel->vertices[i]->coords += avgDispVect * (0.25 - k * 0.025);
		}

		approximateMesh->moveHandles(augMesh, embSkel, optSkel);

		for(int i = 0 ; i < augMesh->N ; i++)
			if(augMesh->freeVertices[i])
				approximateMesh->vertices[i]->coords = augMesh->vertices[i]->coords;
			else
				approximateMesh->vertices[i]->coords = this->vertices[i]->coords;

		approximateMesh->solveWithARAP(augMesh);

		approximateMesh->computeEdgeLengths();
	}
	// ------------------------------------------------------------------------------

	float approximationError = 0.0;

	for(int i = 0 ; i < approximateMesh->N ; i++)
		approximationError += (approximateMesh->vertices[i]->coords - this->vertices[i]->coords).norm();

	printf("Total approximation error: %f\n", approximationError);

	return make_tuple(approximateMesh, optSkel);
}

float Mesh::computeSkelError(Mesh *augMesh, vector<vector<int>> & boneSamples, Skeleton *optSkel)
{
	float totalDist = 0.0;
	vector<tuple<Vector3f, bool>> newSampleRecords = this->computeGtSkelSamples(augMesh);

	for(int i = 0 ; i < boneSamples.size() ; i += 2)
		for(int j = 0 ; j < boneSamples[i].size() ; j++)
		{
			bool isMeaningful = get<1>(newSampleRecords[boneSamples[i][j]]);

			if(isMeaningful)
			{
				Vector3f P = get<0>(newSampleRecords[boneSamples[i][j]]);
				Vector3f U = optSkel->vertices[optSkel->edges[i]->u]->coords;
				Vector3f V = optSkel->vertices[optSkel->edges[i]->v]->coords;

				totalDist += optSkel->computeDistance2LineSeg(U, V, P);
			}
		}

	return totalDist;
}

vector<vector<int>> Mesh::assignSamples2Skel(const Skeleton *embSkel)
{
	vector<vector<int>> boneSamples;

	boneSamples.resize(embSkel->edges.size());

	for(int i = N ; i < vertices.size() ; i++)
	{
		int minIdx;
		float currMin = FLT_MAX;
		Vector3f P = vertices[i]->coords;

		for(int j = 0 ; j < embSkel->edges.size() ; j += 2)
		{
			Vector3f U = embSkel->vertices[embSkel->edges[j]->u]->coords;
			Vector3f V = embSkel->vertices[embSkel->edges[j]->v]->coords;

			float currDist = embSkel->computeDistance2LineSeg(U, V, P);

			if(currDist < currMin)
			{
				currMin = currDist;
				minIdx = j;
			}
		}

		boneSamples[minIdx].push_back(i - N);
		boneSamples[minIdx + 1].push_back(i - N);
	}

	return boneSamples;
}

/*
	neighbouthood[i] contains neighbourhood[i]'th sample (between 0 and allSamples.size()-1). 
*/
vector<int> Mesh::computeSphereNeighbourhood(Vector3f P, float radius)
{
	vector<int> neighbourhood;

	for(int i = N ; i < vertices.size() ; i++)
		if((vertices[i]->coords - P).norm() < radius)
			neighbourhood.push_back(i - N);

	return neighbourhood;
} 

Mesh *Mesh::approximate2gtMeshUnconstrained(Mesh *augMesh)
{
	Mesh *approximateMesh = new Mesh();
	vector<tuple<Vector3f, bool>> newSampleRecords = this->computeGtSkelSamples(augMesh);
	vector<Vector3f> newSampleCoords;

	for(int i = 0 ; i < newSampleRecords.size() ; i++)
		newSampleCoords.push_back(get<0>(newSampleRecords[i]));

	approximateMesh->copyFromMesh(augMesh);

	for(int i = 0 ; i < augMesh->vertices.size() ; i++)
		if(i < augMesh->N && augMesh->freeVertices[i])
			approximateMesh->vertices[i]->coords = augMesh->vertices[i]->coords;
		else if(i < augMesh->N)
			approximateMesh->vertices[i]->coords = this->vertices[i]->coords;
		else
			approximateMesh->vertices[i]->coords = newSampleCoords[i - augMesh->N];

	approximateMesh->solveWithARAP(augMesh);

	approximateMesh->computeEdgeLengths();

	float approximationError = 0.0;

	for(int i = 0 ; i < approximateMesh->N ; i++)
		approximationError += (approximateMesh->vertices[i]->coords - this->vertices[i]->coords).norm();

	printf("Total approximation error: %f\n", approximationError);

	return approximateMesh;
}

vector<tuple<Vector3f, bool>> Mesh::computeGtSkelSamples(Mesh *augMesh)
{
	vector<tuple<Vector3f, bool>> newSampleCoords;

	for(int i = augMesh->N ; i < augMesh->vertices.size() ; i++)
	{
		vector<Vector3f> P;
		vector<Vector3f> Q;

		for(int j = 0 ; j < augMesh->allSupports[i - augMesh->N].size() ; j++)
		{
			P.push_back(augMesh->vertices[augMesh->allSupports[i - augMesh->N][j]]->coords);
			Q.push_back(this->vertices[augMesh->allSupports[i - augMesh->N][j]]->coords);
		}

		/*for(int j = 0 ; j < augMesh->vertices[i]->outgoingEdges.size() ; j++)
		{
			int neighIdx = augMesh->edges[augMesh->vertices[i]->outgoingEdges[j]]->v;

			if(augMesh->edges[neighIdx]->w != 0.0)
			{
				P.push_back(augMesh->vertices[neighIdx]->coords);
				Q.push_back(this->vertices[neighIdx]->coords);
			}
		}*/

		if(P.size() > 0)
		{
			tuple<Matrix3f, Vector3f> transformation = this->computeRigidTransform(P, Q);
			Matrix3f R = get<0>(transformation);
			Vector3f t = get<1>(transformation);
			newSampleCoords.push_back(make_tuple(R * augMesh->vertices[i]->coords + t, true));
		}
		else
			newSampleCoords.push_back(make_tuple(augMesh->vertices[i]->coords, false));
	}

	return newSampleCoords;
}

void Mesh::copyFromMesh(Mesh *mesh)
{
	this->N = mesh->N;
	this->M = mesh->M;

	for(int i = 0 ; i < mesh->vertices.size() ; i++)
	{
		this->vertices.push_back(new Vertex(i, mesh->vertices[i]->coords));
		this->vertices[i]->outgoingEdges = mesh->vertices[i]->outgoingEdges;
	}

	for(int i = 0 ; i < mesh->tris.size() ; i++)
		this->tris.push_back(new Triangle(mesh->tris[i]->id, mesh->tris[i]->v1, mesh->tris[i]->v2, mesh->tris[i]->v3));

	for(int i = 0 ; i < mesh->edges.size() ; i++)
		this->edges.push_back(new Edge(mesh->edges[i]->id, mesh->edges[i]->u, mesh->edges[i]->v, mesh->edges[i]->opp, mesh->edges[i]->faceId));

	this->freeVertices = (bool *) malloc(sizeof(bool) * mesh->vertices.size());
	this->partitions = (int *) malloc(sizeof(int) * mesh->vertices.size());
	this->parametricCoords = (float *) malloc(sizeof(float) * mesh->vertices.size());

	for(int i = 0 ; i < mesh->vertices.size() ; i++)
	{
		this->freeVertices[i] = mesh->freeVertices[i];
		this->partitions[i] = mesh->partitions[i];
		this->parametricCoords[i] = mesh->parametricCoords[i];
	}

	this->ratio = mesh->ratio;
}

void Mesh::getPlotPrimitives(MatrixXd & meshVertices, MatrixXi & meshFaces, MatrixXd & faceColors) const
{
	meshVertices.resize(this->N, 3);
	meshFaces.resize(tris.size(), 3);
	faceColors.resize(tris.size(), 3);

	for(int i = 0 ; i < this->N ; i++)
	{
		meshVertices(i, 0) = vertices[i]->coords(0);
		meshVertices(i, 1) = vertices[i]->coords(1);
		meshVertices(i, 2) = vertices[i]->coords(2);
	}

	for(int i = 0 ; i < tris.size() ; i++)
	{
		meshFaces(i, 0) = tris[i]->v1;
		meshFaces(i, 1) = tris[i]->v2;
		meshFaces(i, 2) = tris[i]->v3;

		if(freeVertices[tris[i]->v1] || freeVertices[tris[i]->v2] || freeVertices[tris[i]->v3])
			faceColors.row(i) << 0.0, 0.7, 0.0;
		else
			faceColors.row(i) << 0.7, 0.7, 0.0;
	}
}

void Mesh::getVerticesFaces(MatrixXd & meshVertices, MatrixXi & meshFaces) const
{
	meshVertices.resize(this->N, 3);
	meshFaces.resize(tris.size(), 3);

	for(int i = 0 ; i < this->N ; i++)
	{
		meshVertices(i, 0) = vertices[i]->coords(0);
		meshVertices(i, 1) = vertices[i]->coords(1);
		meshVertices(i, 2) = vertices[i]->coords(2);
	}

	for(int i = 0 ; i < tris.size() ; i++)
	{
		meshFaces(i, 0) = tris[i]->v1;
		meshFaces(i, 1) = tris[i]->v2;
		meshFaces(i, 2) = tris[i]->v3;
	}
}

void Mesh::getDifferenceColorMap(Mesh *approx, MatrixXd & vertexColors)
{
	vertexColors.resize(this->N, 3);
	vector<float> dispVect;
	vector<float> avgOneRingEdgeLengths;
	float approximationError = 0.0;

	for(int i = 0 ; i < this->N ; i++)
	{
		int numMeshNeighs = 0;
		float disp = (this->vertices[i]->coords - approx->vertices[i]->coords).norm();
		float oneRingAvgEdgeLength = 0.0;

		for(int j = 0 ; j < this->vertices[i]->outgoingEdges.size() ; j++)
		{
			int edgeId = this->vertices[i]->outgoingEdges[j];

			if(edgeId < this->M)
			{
				oneRingAvgEdgeLength += this->edges[edgeId]->length;
				numMeshNeighs++;
			}
		}

		oneRingAvgEdgeLength /= numMeshNeighs;

		dispVect.push_back(disp);
		avgOneRingEdgeLengths.push_back(oneRingAvgEdgeLength);
		approximationError += disp;
	}

	int numFreeVertices = 0;

	for(int i = 0 ; i < approx->N ; i++)
		if(approx->freeVertices[i])
			numFreeVertices++;

	float avgDisp = 0.0;

	for(int i = 0 ; i < approx->N ; i++)
		avgDisp += dispVect[i] / avgOneRingEdgeLengths[i];

	printf("Average displacement of all vertices: %f\n", avgDisp / approx->N);
	printf("Average displacement of only free vertices: %f\n", avgDisp / numFreeVertices);

	for(int i = 0 ; i < this->N ; i++)
	{
		float disp = dispVect[i];
		float red;
		float gre;
		float blu;
		float middle = avgOneRingEdgeLengths[i];

		if(disp < middle)
		{
			red = 0.0;
			gre = disp / middle;
			blu = 1.0 - gre;
		}
		else
		{
			red = disp / middle - 1.0;
			gre = 1.0 - red;
			blu = 0.0;
		}

		red = red > 1.0 ? 1.0 : red;
		gre = gre > 1.0 ? 1.0 : gre;
		blu = blu > 1.0 ? 1.0 : blu;

		vertexColors.row(i) << red, gre, blu;
	}
}

/*
	Returns the number of the initial mesh (non-augmented) vertices. 
*/
unsigned int Mesh::getNumOfMeshVertices(void) const
{
	return this->N;
}

/*
	Returns the number of the initial mesh (non-augmented) edges. 
*/
unsigned int Mesh::getNumOfMeshEdges(void) const
{
	return this->M;
}

/*
	Returns the coordinates of the skeleton vertices (Skeleton samples) in the format libIGL needs. 
*/
void Mesh::getSupportVertices(MatrixXd & supportVertices) const
{
	supportVertices.resize(vertices.size() - this->N, 3);

	for(int i = this->N ; i < vertices.size() ; i++)
	{
		supportVertices(i - this->N, 0) = vertices[i]->coords(0);
		supportVertices(i - this->N, 1) = vertices[i]->coords(1);
		supportVertices(i - this->N, 2) = vertices[i]->coords(2);
	}
}

/*
	Returns the u (dest.) and v (source) of mesh edges. 
*/
vector<tuple<int, int>> Mesh::getMeshEdges(void) const
{
	vector<tuple<int, int>> meshEdges;

	for(int i = 0 ; i < this->M ; i += 2)
		meshEdges.push_back(make_tuple(edges[i]->u, edges[i]->v));

	return meshEdges;
}

/*
	Returns the u (dest.) and v (source) of support edges. 
	u is on the skeleton. 
	v is in the mesh. 
	w is the weight of the edge. 
*/
vector<tuple<int, int, float>> Mesh::getSupportEdges(void) const
{
	vector<tuple<int, int, float>> supportEdges;

	for(int i = this->N ; i < vertices.size() ; i++)
		for(int j = 0 ; j < vertices[i]->outgoingEdges.size() ; j++)
		{
			int meshVert = edges[vertices[i]->outgoingEdges[j]]->v;
			float w = edges[vertices[i]->outgoingEdges[j]]->w;

			supportEdges.push_back(make_tuple(i, meshVert, w));
		}

	return supportEdges;
}

/*
	Returns edge lengths. Length of i'th edge is contained in the i'th position of the returned vector. 
*/
vector<float> Mesh::getEdgeLengths(void) const
{
	vector<float> edgeLengths;

	for(int i = 0 ; i < edges.size() ; i++)
		edgeLengths.push_back(edges[i]->length);

	return edgeLengths;
}

/*
	Returns non-rigidity per 3d cell for each vertex according to the ARAP deformation. 
*/
vector<float> Mesh::getNonRigidity(void) const
{
	return nonRigidityErrors;
}

/*
	Returns weights of the edges. 
*/
vector<float> Mesh::getEdgeWeights(void) const
{
	vector<float> edgeWeights;

	for(int i = 0 ; i < edges.size() ; i++)
		edgeWeights.push_back(edges[i]->w);

	return edgeWeights;
}

/*
	Sets edge weights manually. And then updates vertex weights accordingly. 
*/
void Mesh::setEdgeWeights(vector<float> weights)
{
	for(int i = 0 ; i < edges.size() ; i++)
		edges[i]->w = weights[i];
	
	computeVertexWeights();
}

/*
	Computes mesh volume by projecting each triangle onto xy-plane. 
*/
float Mesh::computeMeshVolume(void)
{
	float volume = 0.0;

	for(int i = 0 ; i < tris.size() ; i++)
	{
		Vector3f P1 = vertices[tris[i]->v1]->coords;
		Vector3f P2 = vertices[tris[i]->v2]->coords;
		Vector3f P3 = vertices[tris[i]->v3]->coords;

		volume += ((P1(2) + P2(2) + P3(2)) / 6.0) * ((P2(0) - P1(0)) * (P3(1) - P1(1)) - (P3(0) - P1(0)) * (P2(1) - P1(1)));
	}

	return volume;
}

/*
	Computes all pairs shortest paths by using Floyd-Warshall algorithm. 
*/
vector<vector<float>> Mesh::computeAPSP(void)
{
	vector<vector<float>> costMatrix(this->N, vector<float>(this->N));

	for(int i = 0 ; i < this->N ; i++)
		for(int j = 0 ; j < this->N ; j++)
			costMatrix[i][j] = FLT_MAX;

	for(int i = 0 ; i < this->M ; i++)
		costMatrix[edges[i]->u][edges[i]->v] = edges[i]->length;

	for(int i = 0 ; i < this->N ; i++)
		costMatrix[i][i] = 0.0;

	for(int k = 0 ; k < this->N ; k++)
		for(int i = 0 ; i < this->N ; i++)
			for(int j = 0 ; j < this->N ; j++)
				if(costMatrix[i][j] > costMatrix[i][k] + costMatrix[k][j])
					costMatrix[i][j] = costMatrix[i][k] + costMatrix[k][j];

	return costMatrix;
}

/*
	All the postprocessing stuff. 
*/
void Mesh::postProcess(void)
{
	int numIter = 10;

	for(int i = 0 ; i < numIter ; i++)
		taubinMeshFairingOneStep();
}

tuple<Matrix3f, Vector3f> Mesh::computeRigidTransform(vector<Vector3f> & P, vector<Vector3f> & Q)
{
	int n = P.size();
	Vector3f P_center(0.0, 0.0, 0.0);
	Vector3f Q_center(0.0, 0.0, 0.0);
	Vector3f t;
	MatrixXf X(3, n);
	MatrixXf Y(3, n);
	MatrixXf S(3, 3);
	MatrixXf R(3, 3);
	MatrixXf W = MatrixXf::Zero(n, n);
	JacobiSVD<MatrixXf> svd;

	for(int i = 0 ; i < n ; i++)
	{
		P_center += P[i];
		Q_center += Q[i];
	}
	P_center /= n;
	Q_center /= n;

	for(int i = 0 ; i < n ; i++)
	{
		X.col(i) = P[i] - P_center;
		Y.col(i) = Q[i] - Q_center;
		W(i, i) = 1.0;
	}

	S = X * W * Y.transpose();
	svd.compute(S, Eigen::ComputeThinU | Eigen::ComputeThinV);
	if(!svd.computeU() || !svd.computeV())
	{
		cout << "Decomposition error." << endl;
		exit(0);
	}

	MatrixXf V = svd.matrixV();
	MatrixXf U = svd.matrixU();
	MatrixXf VUt = V * U.transpose();
	MatrixXf M(3, 3);

	M.setIdentity();

	float VUt_det = VUt.determinant();
	if(VUt_det < 0)
		M(2, 2) = VUt_det;

	R = V * M * U.transpose();
	t = Q_center - R * P_center;

	return make_tuple(R, t);
}

/*
	Samples the skeleton. For each sample we compute closest numOfSupports mesh vertices to it. 
*/
void Mesh::constructAugmentedMeshV1(const Skeleton *skel)
{
	// Take samples on the skeleton with avgEdgeLength between the samples. 
	vector<tuple<Vector3f, Vector3f>> samples = skel->sampleSkeleton(avgEdgeLength);

	// Use samples to augment the mesh. Compute distance to each mesh vertex from each sample. 
	for(int i = 0 ; i < samples.size() ; i++)
	{
		Vector3f sampleCoords = get<0>(samples[i]); 
		vector<tuple<float, int>> distances;

		// Add i'th skeleton vertex to the augmented mesh structure. 
		int supportVertexId = this->N + i;
		vertices.push_back(new Vertex(supportVertexId, sampleCoords));

		for(int j = 0 ; j < this->N ; j++)
		{
			float dist = (vertices[j]->coords - sampleCoords).norm();

			distances.push_back(make_tuple(dist, j));
		}

		sort(distances.begin(), distances.end());

		// Choose closest numOfSupports vertices to i'th sample such that no support vertex is in another's 
		// 1-ring neighbourhood. 
		bool chosen[this->N];	// Chosen mesh vertices as supports. Will be used when not to choose 1-ring neighbourhoods.  
		int numOfSupports = 32;	// TODO: This will be chosen in a smarter way! 
		
		for(int j = 0 ; j < this->N ; j++)
			chosen[j] = false;

		for(int j = 0, k = 0 ; j < numOfSupports ; j++)
		{
			int meshVertexId;

			do
			{
				meshVertexId = get<1>(distances[k]);	// Id of the current surface vertex.
				k++;
			}
			while(chosen[meshVertexId]);
			chosen[meshVertexId] = true;
			 
			for(int m = 0 ; m < vertices[meshVertexId]->outgoingEdges.size() ; m++)
			{
				int dest = edges[vertices[meshVertexId]->outgoingEdges[m]]->v;

				if(dest < this->N)	// Destination should belong to the initial mesh
					chosen[edges[vertices[meshVertexId]->outgoingEdges[m]]->v] = true;
			}

			// Create a support between skeleton vertex i and mesh vertex meshVertexId. 
			int supportEdgeId = this->M + 2 * i * numOfSupports + 2 * j;
			int oppSupportEdgeId = supportEdgeId + 1;

			edges.push_back(new Edge(supportEdgeId, supportVertexId, meshVertexId, oppSupportEdgeId, -1));	// TODO: Should there be a face attached? 
			edges.push_back(new Edge(oppSupportEdgeId, meshVertexId, supportVertexId, supportEdgeId, -1));

			vertices[supportVertexId]->outgoingEdges.push_back(supportEdgeId);
			vertices[meshVertexId]->outgoingEdges.push_back(oppSupportEdgeId);
		}
	}

	computeEdgeLengths();

	computeEdgeWeights();

	// Calculate w's of the support edges. 
	for(int i = this->M ; i < edges.size() ; i += 2)
	{
		edges[i]->w = 1.0 / vertices[edges[i]->u]->outgoingEdges.size();
		edges[i + 1]->w = edges[i]->w;
	}
}

/*
	Samples the skeleton. Adds a support vector between each mesh vertex and closest skeleton sample to it. 
*/
void Mesh::constructAugmentedMeshV2(const Skeleton *skel)
{
	// Take samples on the skeleton with avgEdgeLength between the samples. 
	vector<tuple<Vector3f, Vector3f>> samples = skel->sampleSkeleton(avgEdgeLength);

	for(int i = 0 ; i < samples.size() ; i++)
		vertices.push_back(new Vertex(this->N + i, get<0>(samples[i])));

	for(int i = 0 ; i < this->N ; i++)
	{
		int minIdx = -1;
		float minDist = FLT_MAX;
		Vector3f P_mesh = vertices[i]->coords;

		for(int j = 0 ; j < samples.size() ; j++)
		{
			Vector3f P_skel = get<0>(samples[j]);
			float currDist = (P_mesh - P_skel).norm();

			if(currDist < minDist)
			{
				minDist = currDist;
				minIdx = j;
			}
		}

		// Create a support between mesh vertex and skeleton vertex. 
		int supportVertexId = this->N + minIdx;
		int supportEdgeId = this->M + 2 * i;
		int oppSupportEdgeId = supportEdgeId + 1;

		edges.push_back(new Edge(supportEdgeId, supportVertexId, i, oppSupportEdgeId, -1));	// TODO: Should there be a face attached? 
		edges.push_back(new Edge(oppSupportEdgeId, i, supportVertexId, supportEdgeId, -1));

		vertices[supportVertexId]->outgoingEdges.push_back(supportEdgeId);
		vertices[i]->outgoingEdges.push_back(oppSupportEdgeId);
	}

	computeEdgeLengths();

	computeEdgeWeights();

	// Calculate w's of the support edges. 
	for(int i = this->M ; i < edges.size() ; i += 2)
	{
		edges[i]->w = 1.0 / vertices[edges[i]->u]->outgoingEdges.size();
		edges[i + 1]->w = edges[i]->w;
	}
}

void Mesh::constructAugmentedMeshContour(const Skeleton *skel)
{
	// Take samples on the skeleton with avgEdgeLength between the samples. 
	printf("Augmented mesh is being constructed...\n");
	printf("\tSampling skeleton...\n");
	vector<tuple<Vector3f, Vector3f>> samples = skel->sampleSkeleton(avgEdgeLength);
	printf("\tSkeleton sampled.\n");

	auto start = std::chrono::steady_clock::now();

	printf("\tSupport vertices and edges is being computed...\n");
	for(int i = 0 ; i < samples.size() ; i++)
	{
		vector<int> supports;	// Support edges emanating from i'th skeleton vertex. 
		Vector3f P_p = get<0>(samples[i]);	// Point on the plane. A skeleton vertex. 
		Vector3f n_p = get<1>(samples[i]);	// Plane normal. 

		for(int j = 0 ; j < this->N ; j++)
		{
			Vector3f P_0 = vertices[j]->coords;	// Mesh vertex coordinate. 
			float d = fabs((P_0 - P_p).dot(n_p));

			if(d <= avgEdgeLength / 2.0)	// Candidate, further checks must be done. 
			{
				Vector3f dir = (P_0 - P_p).normalized();	// Ray direction. 
				Vector3f org = P_p + 0.0001 * dir;	// Ray origin. 

				Ray skeletonRay(org, dir);
				tuple<Vector3f, float> retVal = bvh->intersect(skeletonRay, bvh);
				Vector3f intPoint = get<0>(retVal);

				if(intPoint(0) != FLT_MAX || intPoint(1) != FLT_MAX || intPoint(2) != FLT_MAX)
				{
					float d_int = (intPoint - org).norm();	// Distance between mesh intersection point and ray origin. 
					float l_ray = (P_0 - org).norm();	// Expected distance. 

					if(d_int + 0.0001 < l_ray)	// There is an 'obstacle', reject. 
						continue;
				}

				supports.push_back(j);
			}
		}

		allSupports.push_back(supports);

		int supportVertIdx = this->N + i;

		vertices.push_back(new Vertex(supportVertIdx, P_p));
	}

	auto end = std::chrono::steady_clock::now();

	std::chrono::duration<double> elapsedSeconds = end - start;

	printf("\tSupport vertices and edges computed. (%lf seconds)\n", elapsedSeconds.count());

	freeVertices = (bool *) malloc(sizeof(bool) * vertices.size());
	partitions = (int *) malloc(sizeof(int) * vertices.size());
	parametricCoords = (float *) malloc(sizeof(float) * vertices.size());

	correspondMesh(skel);

	markFreeVertices(ratio, freeVertices, skel);

	// Paint regions. 

	start = std::chrono::steady_clock::now();

	printf("\tValid support edges is being computed...\n");
	printf("\t\tPainting regions...\n");

	int currColor = 1;

	for(int i = 0 ; i < vertices.size() ; i++)
		if(vertices[i]->region < 0)
		{
			queue<int> bfsQue;	// Standard FIFO queue. Implements Breadth-First-Search (BFS). 

			vertices[i]->region = currColor;

			bfsQue.push(i);

			while(!bfsQue.empty())
			{
				int currVertex = bfsQue.front();

				bfsQue.pop();

				for(int j = 0 ; j < vertices[currVertex]->outgoingEdges.size() ; j++)
				{
					int dest = edges[vertices[currVertex]->outgoingEdges[j]]->v;

					if(vertices[dest]->region < 0)
					{
						vertices[dest]->region = currColor;

						bfsQue.push(dest);
					}
				}
			}

			currColor++;
		}
	printf("\t\tRegions painted.\n");

	int supportEdgeIdx = this->M;

	for(int i = 0 ; i < allSupports.size() ; i++)
	{	
		int interestedRegion;
		int supportVertIdx = this->N + i;
		float currMin = FLT_MAX;
		Vector3f P_p = vertices[supportVertIdx]->coords;

		for(int j = 0 ; j < allSupports[i].size() ; j++)
		{
			float currDist = (vertices[allSupports[i][j]]->coords - P_p).norm();

			if(currDist < currMin && vertices[allSupports[i][j]]->region != 0)
			{
				currMin = currDist;
				interestedRegion = vertices[allSupports[i][j]]->region;
			}
		}

		for(int j = 0 ; j < allSupports[i].size() ; j++)
			if(vertices[allSupports[i][j]]->region == interestedRegion)
			{
				int oppSupportEdgeIdx = supportEdgeIdx + 1;

				edges.push_back(new Edge(supportEdgeIdx, supportVertIdx, allSupports[i][j], oppSupportEdgeIdx, -1));	// TODO: Should there be a face attached? 
				edges.push_back(new Edge(oppSupportEdgeIdx, allSupports[i][j], supportVertIdx, supportEdgeIdx, -1));

				vertices[supportVertIdx]->outgoingEdges.push_back(supportEdgeIdx);
				vertices[allSupports[i][j]]->outgoingEdges.push_back(oppSupportEdgeIdx);

				supportEdgeIdx += 2;
			}
	}

	end = std::chrono::steady_clock::now();

	elapsedSeconds = end - start;

	printf("\tValid support edges computed. (%lf seconds)\n", elapsedSeconds.count());

	computeEdgeLengths();

	printf("\tWeights are being computed.\n");
	computeEdgeWeights();

	// Calculate w's of the support edges. 
	for(int i = 0 ; i < this->N ; i++)
	{
		int numOfMeshEdges = 0;
		float tempSum = 0.0;
		float w_is = 0.0;
		float minDist = FLT_MAX;
		int closestInd = -1;

		for(int j = 0 ; j < vertices[i]->outgoingEdges.size() ; j++)
		{
			int edgeIdx = vertices[i]->outgoingEdges[j];

			if(edgeIdx < this->M)
			{
				w_is += edges[edgeIdx]->w;
				numOfMeshEdges++;
			}
			else
			{
				float currDist = edges[edgeIdx]->length;

				if(currDist < minDist)
				{
					minDist = currDist;
					closestInd = edgeIdx;
				}
			}
		}
		
		w_is /= numOfMeshEdges;

		if(closestInd != -1)
		{
			edges[closestInd]->w = w_is;
			edges[closestInd - 1]->w = w_is;
		}
	}
	
	int numAdditionalVerts = 0;
	for(int i = this->N ; i < vertices.size() ; i++)
	{
		for(int j = 0 ; j < vertices[i]->outgoingEdges.size() ; j++)
		{
			int edgeIdx = vertices[i]->outgoingEdges[j];

			if(freeVertices[edges[edgeIdx]->v] && edges[edgeIdx]->w != 0.0)
			{
				numAdditionalVerts++;
				break;
			}
		}
	}

	printf("\tWeights computed.\n");
	printf("Augmented mesh constructed (%d support vertices).\n", numAdditionalVerts);
}

// Standard arap implementation does not need any additional structure. 
void Mesh::constructMeshForARAP(const Skeleton *skel)
{
	freeVertices = (bool *) malloc(sizeof(bool) * vertices.size());
	partitions = (int *) malloc(sizeof(int) * vertices.size());
	parametricCoords = (float *) malloc(sizeof(float) * vertices.size());

	correspondMesh(skel);

	markFreeVertices(ratio, freeVertices, skel);

	computeEdgeLengths();

	computeEdgeWeights();
}

/*
	Augments the mesh by using Zhang's approach. 
*/
void Mesh::constructAugmentedMeshZhang(const Skeleton *skel)
{
	int raysPerSample = 6;
	bool duplicateCheck[this->N];
	vector<tuple<Vector3f, vector<int>>> meshUpdateStructure; 
	vector<tuple<Vector3f, vector<Vector3f>>> orthRays = skel->computeRays(avgEdgeLength, raysPerSample);

	for(int i = 0 ; i < orthRays.size() ; i++)
	{
		vector<int> meshVerts;	// Nearest vertices to a skeleton vertex on the mesh surface. 
		Vector3f org = get<0>(orthRays[i]);	// Ray origin
		vector<Vector3f> dirs = get<1>(orthRays[i]);	// Ray directions 

		for(int j = 0 ; j < raysPerSample ; j++)
		{
			Ray skeletonRay(org, dirs[j]);
			Vector3f intPoint = rayMeshIntersection(skeletonRay);

			if(intPoint(0) != 0.0 || intPoint(1) != 0.0 || intPoint(2) != 0.0)
			{
				int nearestIdx = nearestVertex(intPoint);

				if(!duplicateCheck[nearestIdx])
				{
					meshVerts.push_back(nearestIdx);
					duplicateCheck[nearestIdx] = true;
				}
			}
		}

		for(int j = 0 ; j < meshVerts.size() ; j++)
			duplicateCheck[meshVerts[j]] = false;

		if(meshVerts.size() > 2)
			meshUpdateStructure.push_back(make_tuple(org, meshVerts));
	}

	int edgeIdx = this->M;
	for(int i = 0 ; i < meshUpdateStructure.size() ; i++)
	{
		int vertIdx = this->N + i;
		Vector3f coords = get<0>(meshUpdateStructure[i]);
		vector<int> neighbours = get<1>(meshUpdateStructure[i]);

		vertices.push_back(new Vertex(vertIdx, coords));

		for(int j = 0 ; j < neighbours.size() ; j++)
		{
			edges.push_back(new Edge(edgeIdx, vertIdx, neighbours[j], edgeIdx + 1, -1));
			edges.push_back(new Edge(edgeIdx + 1, neighbours[j], vertIdx, edgeIdx, -1));

			vertices[vertIdx]->outgoingEdges.push_back(edgeIdx);
			vertices[neighbours[j]]->outgoingEdges.push_back(edgeIdx + 1);

			edgeIdx += 2;
		}
	}

	freeVertices = (bool *) malloc(sizeof(bool) * vertices.size());
	partitions = (int *) malloc(sizeof(int) * vertices.size());
	parametricCoords = (float *) malloc(sizeof(float) * vertices.size());

	correspondMesh(skel);

	markFreeVertices(ratio, freeVertices, skel);

	computeEdgeLengths();

	computeEdgeWeights();

	// Calculate w's of the support edges. 
	for(int i = this->M ; i < edges.size() ; i += 2)
	{
		//edges[i]->w = 1.0 / vertices[edges[i]->u]->outgoingEdges.size();
		edges[i]->w = 1.0;
		edges[i + 1]->w = edges[i]->w;
	}

	for(int i = this->N ; i < vertices.size() ; i++)
	{
		unsigned int numOutEdges = vertices[i]->outgoingEdges.size();
		// Coordinates of the skeleton vertex. 
		Vector3f V_s = vertices[i]->coords;

		for(int j = 0 ; i < numOutEdges ; j++)
		{
			int edgeIdx = edges[vertices[i]->outgoingEdges[j]]->v;

			// Coordinates of the "current" vertex. 
			Vector3f V_c = vertices[edgeIdx]->coords;

			// Coordinates of the "previous" vertex. 
			Vector3f V_p = vertices[edges[vertices[i]->outgoingEdges[(j + (numOutEdges - 1)) % numOutEdges]]->v]->coords;

			// Coordinates of the "next" vertex. 
			Vector3f V_n = vertices[edges[vertices[i]->outgoingEdges[(j + 1) % numOutEdges]]->v]->coords;

			// Compute alpha. 
			float cosalpha = ((V_c - V_p).normalized()).dot((V_s - V_p).normalized());

			if(cosalpha > 1.0)
				cosalpha = 1.0;
			if(cosalpha < -1.0)
				cosalpha = -1.0;

			float alpha = acos(cosalpha);

			// Compute beta. 
			float cosbeta = ((V_c - V_n).normalized()).dot((V_s - V_n).normalized());

			if(cosbeta > 1.0)
				cosbeta = 1.0;
			if(cosbeta < -1.0)
				cosbeta = -1.0;

			float beta = acos(cosbeta);

			edges[edgeIdx]->w = 0.5 * (1.0 / tan(alpha) + 1.0 / tan(beta));
			edges[edgeIdx + 1]->w = edges[edgeIdx]->w;
		}
	}
}

void Mesh::moveHandles(Mesh *sMesh, const Skeleton *sSkel, const Skeleton *tSkel)
{
	// Index of outer vector is the edge id. 
	// Inner vectors hold vertex id's corresponded with the edge id which is outer vector index. 
	vector<vector<int>> partitionsByEdges(sSkel->edges.size());

	for(int i = 0 ; i < vertices.size() ; i++)
	{
		partitionsByEdges[partitions[i]].push_back(i);
		partitionsByEdges[partitions[i] + 1].push_back(i);
	}

	int trunkIdx = sSkel->findTrunk();	// Find the trunk edge id (The edge whose deg(u)+deg(v) is maximum).
	float optRAAAngles[sSkel->edges.size()];	// Optimum Rotation-Around Axis (RAA) Angles. 
	Vector3f handles[sSkel->edges.size()];	// Handles as unit vectors. The end of a handle is always at the half-edge. 
	Vector3f transferredHandles[sSkel->edges.size()];	//  Handles as unit vectors transferred into the target skeleton. 

	initializeHandles(sSkel, handles);	// Compute handles for the source skeleton. 
	optRAAAngles[trunkIdx] = computeOptRAAForTrunk(trunkIdx, sSkel, tSkel);	// Computation of optimum RAA for trunk is different than other RAAs. 
	optRAAAngles[sSkel->edges[trunkIdx]->opp] = optRAAAngles[trunkIdx];

	bool visitedEdges[sSkel->edges.size()];	// To keep track of visited edges in queue. 
	float axisAlignedAngles[sSkel->edges.size()];	// Axis-Aligned Angles relative to the parent.

	for(int i = 0 ; i < sSkel->edges.size() ; i++)
		visitedEdges[i] = false;
	visitedEdges[trunkIdx] = true;	// Trunk edge has been processed. 
	visitedEdges[sSkel->edges[trunkIdx]->opp] = true;

	queue<int> bfsQue;	// Standard FIFO queue. Implements Breadth-First-Search (BFS). 

	bfsQue.push(trunkIdx);	// Initialize queue with trunk edge. 
	while(!bfsQue.empty())
	{
		int parentEdge = bfsQue.front();
		vector<int> childrenEdges;

		// Transfer the partition corresponding to parentEdge (translation+edgeAlignment+rotationAlignment). 
		for(int i = 0 ; i < partitionsByEdges[parentEdge].size() ; i++)
		{
			int vertIdx = partitionsByEdges[parentEdge][i];
			Vector3f vertCoords = sMesh->vertices[vertIdx]->coords;

			vertCoords = applyEdgeTranslation(vertCoords, parentEdge, sSkel, tSkel);
			vertCoords = applyEdgeAlignment(vertCoords, parentEdge, sSkel, tSkel);
			//vertCoords = rotateAroundEdge(vertCoords, optRAAAngles[parentEdge], parentEdge, tSkel);

			this->vertices[vertIdx]->coords = vertCoords;
		}

		// Do the same for handles[parentEdge] and set targetHandles[parentEdge].
		Vector3f Up = sSkel->vertices[sSkel->edges[parentEdge]->u]->coords;
		Vector3f Vp = sSkel->vertices[sSkel->edges[parentEdge]->v]->coords;
		Vector3f Hp = (Up + Vp) / 2 + handles[parentEdge];

		Hp = applyEdgeTranslation(Hp, parentEdge, sSkel, tSkel);
		Hp = applyEdgeAlignment(Hp, parentEdge, sSkel, tSkel);
		Hp = rotateAroundEdge(Hp, optRAAAngles[parentEdge], parentEdge, tSkel);
		Up = tSkel->vertices[tSkel->edges[parentEdge]->u]->coords;
		Vp = tSkel->vertices[tSkel->edges[parentEdge]->v]->coords;
		transferredHandles[parentEdge] = Hp - (Up + Vp) / 2;
		transferredHandles[tSkel->edges[parentEdge]->opp] = transferredHandles[parentEdge];

		bfsQue.pop();	// Parent edge has been processed. 

		// Arrange children edges. 
		int u = sSkel->edges[parentEdge]->u;
		int v = sSkel->edges[parentEdge]->v;

		for(int i = 0 ; i < sSkel->vertices[u]->outgoingEdges.size() ; i++)
		{
			int childEdge = sSkel->vertices[u]->outgoingEdges[i];

			if(!visitedEdges[childEdge])
				childrenEdges.push_back(childEdge);
		}
		for(int i = 0 ; i < sSkel->vertices[v]->outgoingEdges.size() ; i++)
		{
			int childEdge = sSkel->vertices[v]->outgoingEdges[i];

			if(!visitedEdges[childEdge])
				childrenEdges.push_back(childEdge);
		}

		for(int i = 0 ; i < childrenEdges.size() ; i++)
		{
			int childEdge = childrenEdges[i];
		
			bfsQue.push(childEdge);
			visitedEdges[childEdge] = true;
			visitedEdges[sSkel->edges[childEdge]->opp] = true;

			// Compute Axis-Aligned Angle (AAA) between parentEdge and childEdge in the source skeleton.
			// This is the AAA that should be between the parent edge and child edge in the target skeleton. 
			//float AAAngleSource = computeAxisAlignedAngle(childEdge, parentEdge, sSkel, handles); 

			// Transfer handle[childEdge] (translation+edgeAlignment). 
			// What is the AAA between parent and child when we move the handle to the target with only translation+edgeAlignment? 
			Vector3f Uc = sSkel->vertices[sSkel->edges[childEdge]->u]->coords;
			Vector3f Vc = sSkel->vertices[sSkel->edges[childEdge]->v]->coords;
			Vector3f Hc = (Uc + Vc) / 2 + handles[childEdge];

			Hc = applyEdgeTranslation(Hc, childEdge, sSkel, tSkel);
			Hc = applyEdgeAlignment(Hc, childEdge, sSkel, tSkel);
			Uc = tSkel->vertices[tSkel->edges[childEdge]->u]->coords;
			Vc = tSkel->vertices[tSkel->edges[childEdge]->v]->coords;
			transferredHandles[childEdge] = Hc - (Uc + Vc) / 2;
			transferredHandles[tSkel->edges[childEdge]->opp] = transferredHandles[childEdge];

			// Answer to the question above.
			//float AAAngleTarget = computeAxisAlignedAngle(childEdge, parentEdge, tSkel, transferredHandles); 

			// Compute misalignment between parentEdge and transferred child edge in the target skeleton. 
			//float misalignment = AAAngleSource - AAAngleTarget;
			//float misalignment = M_PI;

			// Set optRAAAngles[childEdge] to misalignment. 
			//childEdge = 7;
			//optRAAAngles[childEdge] = misalignment;
			//optRAAAngles[tSkel->edges[childEdge]->opp] = misalignment;
			//printf("%d - %d - %f - %f - %f\n", childEdge/2, parentEdge/2, AAAngleSource, AAAngleTarget, misalignment);
		}
	}
}

void Mesh::solveWithARAP(Mesh *sMesh)
{
	vector<Matrix3f> R(vertices.size());	// Rotations (for each cell). 
	vector<Triplet<float>> triplets;
	SparseMatrix<float> W(vertices.size(), vertices.size());	// Cotangent weights. 
	SparseMatrix<float> L(vertices.size(), vertices.size());	// Discrete Laplace-Beltrami operator. 

	// Compute w_ij's and L.
	for(int i = 0 ; i < sMesh->vertices.size() ; i++)
	{
		for(int j = 0 ; j < sMesh->vertices[i]->outgoingEdges.size() ; j++)
		{
			int neighVert = sMesh->edges[sMesh->vertices[i]->outgoingEdges[j]]->v;
			float w_ij = sMesh->edges[sMesh->vertices[i]->outgoingEdges[j]]->w;

			W.insert(i, neighVert) = w_ij;

			if(freeVertices[i])	// Check if i is not handle.
				triplets.push_back(Triplet<float>(i, neighVert, -w_ij));
		}

		if(freeVertices[i])	// Check if i is not handle. 
			triplets.push_back(Triplet<float>(i, i, sMesh->vertices[i]->w));	// i is not handle. 
		else
			triplets.push_back(Triplet<float>(i, i, 1.0));	// i is handle. 
	}

	// Arrange L. 
	L.setFromTriplets(triplets.begin(), triplets.end());

	auto startFac = std::chrono::steady_clock::now();

	printf("System matrix factorization stage...\n");

	// Alternating optimization. 
	//BiCGSTAB<SparseMatrix<float>> solver(L);	// BiConjugate Gradients STABilized Method 
	SparseLU<SparseMatrix<float>, COLAMDOrdering<int>> solverLU;

	solverLU.analyzePattern(L);
	solverLU.factorize(L);

	auto endFac = std::chrono::steady_clock::now();

	std::chrono::duration<double> elapsedSecondsFac = endFac - startFac;

	printf("System matrix factorized by LU decomposition (%lf seconds). \n", elapsedSecondsFac.count());

	int currStep = 1;
	int numThreads = thread::hardware_concurrency();
	vector<int> vertsPerThread(numThreads, sMesh->vertices.size() / numThreads);
	vector<tuple<int, int>> threadIndices;
	vector<thread> threads(numThreads);

	for(int residual = sMesh->vertices.size() % numThreads, i = 0 ; residual > 0 ; residual--, i++)
		vertsPerThread[i]++;

	printf("Number of threads available: %d\n", numThreads);

	int head = 0;
	int tail = vertsPerThread[0];

	threadIndices.push_back(make_tuple(head, tail));

	for(int i = 0 ; i < numThreads - 1 ; i++)
	{
		head = tail;
		tail = head + vertsPerThread[i + 1];

		threadIndices.push_back(make_tuple(head, tail));
	}

	vector<vector<int>> neighbourhoodInfo(sMesh->vertices.size());

	for(int i = 0 ; i < sMesh->vertices.size() ; i++)
		for(int j = 0 ; j < sMesh->vertices[i]->outgoingEdges.size() ; j++)
			if(sMesh->edges[sMesh->vertices[i]->outgoingEdges[j]]->v < sMesh->N)
				neighbourhoodInfo[i].push_back(sMesh->edges[sMesh->vertices[i]->outgoingEdges[j]]->v);
		else
			neighbourhoodInfo[i].push_back(sMesh->edges[sMesh->vertices[i]->outgoingEdges[j]]->v);

	printf("Optimization stage...\n");

	auto start = std::chrono::steady_clock::now();

	while(true)
	{
		float totalError;
		float previousTotalError;
		
#ifndef SR_ARAP
		// Compute R_i's for ARAP and SD-ARAP. 
		vector<thread> threads(numThreads);

		for(int i = 0 ; i < numThreads ; i++)
			threads[i] = thread(&Mesh::computeOptimumRotation_, this, std::cref(sMesh), std::ref(W), std::ref(neighbourhoodInfo), std::ref(R), get<0>(threadIndices[i]), get<1>(threadIndices[i]));

		for(int i = 0 ; i < numThreads ; i++)
			threads[i].join();
#endif

#ifdef SR_ARAP
		// Compute R_i's for SR-ARAP. 
		for(int i = 0 ; i < sMesh->vertices.size() ; i++)
			R[i] = computeOptimumSmoothRotation(sMesh, i, W, R);
#endif

		// Compute p_i_'s.
		updatePositions(sMesh, R, solverLU, W, L, freeVertices);

		previousTotalError = totalError;

		totalError = computeError(sMesh, R, W, false);

		//if(currStep > 1 && (previousTotalError - totalError) / previousTotalError < 0.00001)
		//	break;

		currStep++;

		if(currStep > 100) break;
	}

	auto end = std::chrono::steady_clock::now();

	std::chrono::duration<double> elapsedSeconds = end - start;

	printf("Optimization stage ended. (%d steps) (%lf seconds)\n", currStep, elapsedSeconds.count());

	/*for(int i = 0 ; i < this->N ; i++)
	{
		Matrix3f R_mesh;

		computeOptimumRotation(sMesh, i, P, W, true, R_mesh);

		float errorPerCell = computePerCellError(sMesh, i, R_mesh, W, true);
		
		this->nonRigidityErrors.push_back(errorPerCell);
	}*/
}

void Mesh::computeOptimumRotation_(const Mesh *sMesh, SparseMatrix<float> & W, vector<vector<int>> & neighbourhoodInfo, vector<Matrix3f> & R, int b, int e) const
{
	for(int i = b ; i < e ; i++)
		computeOptimumRotation(sMesh, i, W, neighbourhoodInfo, R[i]);
}

/*
	Computes weights of the edges only for the mesh. 
*/
void Mesh::computeEdgeWeights(void)
{
	// Compute weights of the edges. 
	for(int i = 0 ; i < this->M ; i++)
	{
		if(edges[i]->w != 0.0)
			continue;

		float cot_alpha = computeCotan(i);
		float cot_beta = computeCotan(edges[i]->opp);

		edges[i]->w = 0.5 * (cot_alpha + cot_beta);
		edges[edges[i]->opp]->w = edges[i]->w;
	}
}

/*
	Computes weights of the vertices. 
*/
void Mesh::computeVertexWeights(void)
{
	for(int i = 0 ; i < vertices.size() ; i++)
	{
		vertices[i]->w = 0.0;

		for(int j = 0 ; j < vertices[i]->outgoingEdges.size() ; j++)
			vertices[i]->w += edges[vertices[i]->outgoingEdges[j]]->w;
	}
}

Vector3f Mesh::rayMeshIntersection(Ray & ray) const
{
	float minDist = FLT_MAX;
	Vector3f closestPoint(FLT_MAX, FLT_MAX, FLT_MAX);
	tuple<Vector3f, float> retVal;

	for(int i = 0 ; i < tris.size() ; i++)
	{
		retVal = tris[i]->rayTriangleIntersection(ray);

		float t = get<1>(retVal);

		if(t == FLT_MAX)
			continue;

		Vector3f currPoint = get<0>(retVal);
		float currDist = (currPoint - ray.o).norm();

		if(currDist < minDist)
		{
			minDist = currDist;
			closestPoint = currPoint;
		}
	}

	return closestPoint;
}

/*
	Returns the index of the mesh vertex closest to P. 
*/
int Mesh::nearestVertex(Vector3f P) const
{
	int minIdx = -1;
	float minDist = FLT_MAX;

	for(int i = 0 ; i < vertices.size() ; i++)
	{
		float currDist = (vertices[i]->coords - P).norm();

		if(currDist < minDist)
		{
			minDist = currDist;
			minIdx = i;
		}
	}

	return minIdx;
}

/*
	Updates edge lengths. 
*/
void Mesh::computeEdgeLengths(void)
{
	for(int i = 0 ; i < edges.size() ; i++)
	{
		int u = edges[i]->u;
		int v = edges[i]->v;

		edges[i]->length = (vertices[u]->coords - vertices[v]->coords).norm();
	}
}

/*
	Sets normals at the vertices of the mesh as unit vector. 
	Angles between incident edges are used as weights. 

*/
void Mesh::setVerticeNormals(void)
{
	for(int i = 0 ; i < this->N ; i++)
	{
		Vector3f sumOfNormals(0.0, 0.0, 0.0);

		for(int j = 0 ; j < vertices[i]->outgoingEdges.size() ; j++)
		{
			int triId = edges[vertices[i]->outgoingEdges[j]]->faceId;
			Vector3f triNormal = tris[triId]->normal;

			int UIdx = edges[vertices[i]->outgoingEdges[j]]->u;
			int V1Idx = edges[vertices[i]->outgoingEdges[j]]->v;
			int V2Idx = tris[triId]->getThirdVertex(UIdx, V1Idx);
			Vector3f U = vertices[UIdx]->coords;
			Vector3f V1 = vertices[V1Idx]->coords;
			Vector3f V2 = vertices[V2Idx]->coords;
			Vector3f e1 = (V1 - U).normalized();
			Vector3f e2 = (V2 - U).normalized();

			sumOfNormals += acos(e1.dot(e2)) * triNormal;
		}

		vertices[i]->normal = sumOfNormals.normalized();
	}
}

/*
	Computes axis-aligned-angle between two edges, which have a parent-child relationship. 
*/
float Mesh::computeAxisAlignedAngle(int childEdge, int parentEdge, const Skeleton *skel, Vector3f *handles) const
{
	int u_c = skel->edges[childEdge]->u;
	int v_c = skel->edges[childEdge]->v;
	int u_p;
	int v_p;

	if(u_c == skel->edges[parentEdge]->v)
	{
		u_p = skel->edges[parentEdge]->u;
		v_p = skel->edges[parentEdge]->v;
	}
	else
	{
		u_p = skel->edges[parentEdge]->v;
		v_p = skel->edges[parentEdge]->u;
	}

	Vector3f U_c = skel->vertices[u_c]->coords;
	Vector3f V_c = skel->vertices[v_c]->coords;
	Vector3f U_p = skel->vertices[u_p]->coords;
	Vector3f V_p = skel->vertices[v_p]->coords;
	Vector3f wp = (V_p - U_p).normalized();
	Vector3f wc = (V_c - U_c).normalized();

	Vector3f pivot = (wc.cross(wp)).normalized();
	float theta = acos(wc.dot(wp));
	Vector3f hc = rotate(handles[childEdge], theta, pivot);
	Vector3f hp = handles[parentEdge];

	return acos(hc.dot(hp));
}

/*
	Computes handle coordinates. The end of a handle is always at the half-edge. 
*/
void Mesh::initializeHandles(const Skeleton *sSkel, Vector3f *handles) const
{
	Vector3f P(10.0, 10.0, 10.0);	// Dummy point. 

	for(int i = 0 ; i < sSkel->edges.size() ; i += 2)	// Only even edges, the others are just opposites.
	{
		int u = sSkel->edges[i]->u; 
		int v = sSkel->edges[i]->v;
		Vector3f U = sSkel->vertices[u]->coords;
		Vector3f V = sSkel->vertices[v]->coords;
		Vector3f edgeDirection = (V - U).normalized();
		Vector3f P_ = U + edgeDirection * (P - U).dot(edgeDirection);
		Vector3f n = (P - P_).normalized();

		handles[i] = n;
		handles[i + 1] = n;
	} 
}

/*
	Computes optimum RAA for trunk. Return value is in terms of radians around CCW direciton. 
	It will be used as the seed for RAAs of other parts of the mesh. 
*/
float Mesh::computeOptRAAForTrunk(int trunkIdx, const Skeleton *sSkel, const Skeleton *tSkel) const
{
	int u = sSkel->edges[trunkIdx]->u;
	int v = sSkel->edges[trunkIdx]->v;
	vector<int> neighs;
	vector<Vector3f> P_;
	float optRotation;	// CCW. 

	for(int i = 0 ; i < sSkel->vertices[u]->outgoingEdges.size() ; i++)
	{
		int uv = sSkel->vertices[u]->outgoingEdges[i];

		if(sSkel->edges[uv]->v != v)
			neighs.push_back(uv);
	}

	for(int i = 0 ; i < sSkel->vertices[v]->outgoingEdges.size() ; i++)
	{
		int vv = sSkel->vertices[v]->outgoingEdges[i];

		if(sSkel->edges[vv]->v != u)
			neighs.push_back(vv);
	}

	for(int i = 0 ; i < neighs.size() ; i++)
	{
		int vertId = sSkel->edges[neighs[i]]->v;
		Vector3f P = sSkel->vertices[vertId]->coords;

		P = applyEdgeTranslation(P, trunkIdx, sSkel, tSkel);
		P = applyEdgeAlignment(P, trunkIdx, sSkel, tSkel);

		P_.push_back(P);
	}

	float maxScore = 0.0;
	float eps = 0.000005;

	for(int i = 0 ; i < 360 ; i++)
	{
		float currScore = computeE_phi(P_, neighs, tSkel);

		if(currScore > maxScore + eps)
		{
			maxScore = currScore;
			optRotation = ((float) i) / 360.0 * 2 * M_PI;
		}

		for(int j = 0 ; j < P_.size() ; j++)
			P_[j] = rotateAroundEdge(P_[j], 2.0 * M_PI / 360.0, trunkIdx, tSkel);
	}

	return optRotation;
}

float Mesh::computeE_phi(vector<Vector3f> setSrc, vector<int> neighs, const Skeleton *tSkel) const
{
	float score = 0.0;
	vector<Vector3f> e_s;
	vector<Vector3f> e_t;

	for(int i = 0 ; i < neighs.size() ; i++)
	{
		int u = tSkel->edges[neighs[i]]->u;
		int v = tSkel->edges[neighs[i]]->v;

		e_s.push_back((setSrc[i] - tSkel->vertices[u]->coords).normalized());
		e_t.push_back((tSkel->vertices[v]->coords - tSkel->vertices[u]->coords).normalized());
	}

	for(int i = 0 ; i < e_s.size() ; i++)
		score += e_s[i].dot(e_t[i]);

	return score;
}

/*
	Rotates P around edge with edgeId in the target skeleton with angle theta in CCW direction. 
*/
Vector3f Mesh::rotateAroundEdge(Vector3f P, float theta, int edgeId, const Skeleton *tSkel) const
{
	Vector3f U = tSkel->vertices[tSkel->edges[edgeId]->u]->coords;
	Vector3f V = tSkel->vertices[tSkel->edges[edgeId]->v]->coords;
	Vector3f r = (V - U).normalized();

	P -= U;
	P = rotate(P, theta, r);
	P += U;

	return P;
}

/*
	Rotates point P around r in CCW direction with theta (in radians) angle. 
*/
Vector3f Mesh::rotate(Vector3f P, float theta, Vector3f r) const
{
	Vector3f Q(0.0, 0.0, 0.0);

	float costheta = cos(theta);
	float sintheta = sqrt(1.0 - costheta * costheta);

	Q(0) += (costheta + (1 - costheta) * r(0) * r(0)) * P(0);
   	Q(0) += ((1 - costheta) * r(0) * r(1) - r(2) * sintheta) * P(1);
	Q(0) += ((1 - costheta) * r(0) * r(2) + r(1) * sintheta) * P(2);

	Q(1) += ((1 - costheta) * r(0) * r(1) + r(2) * sintheta) * P(0);
	Q(1) += (costheta + (1 - costheta) * r(1) * r(1)) * P(1);
	Q(1) += ((1 - costheta) * r(1) * r(2) - r(0) * sintheta) * P(2);

	Q(2) += ((1 - costheta) * r(0) * r(2) - r(1) * sintheta) * P(0);
	Q(2) += ((1 - costheta) * r(1) * r(2) + r(0) * sintheta) * P(1);
	Q(2) += (costheta + (1 - costheta) * r(2) * r(2)) * P(2);

	return Q;
}

/*
	Save the mesh as .off file. 
*/
void Mesh::saveAsOFF(const char *filePath) const
{
	FILE *pFile = fopen(filePath, "w");

	fprintf(pFile, "OFF\n");
	fprintf(pFile, "%lu %lu 0\n", this->N, tris.size());

	for(int i = 0 ; i < this->N ; i++)
		fprintf(pFile, "%f %f %f\n", vertices[i]->coords(0), vertices[i]->coords(1), vertices[i]->coords(2));

	for(int i = 0 ; i < tris.size() ; i++)
		fprintf(pFile, "3 %d %d %d\n", tris[i]->v1, tris[i]->v2, tris[i]->v3);

	fclose(pFile);
}

/*
	Saves the partition induced by edge with edgeId as OFF file. 
	TODO: Implement triangles also, current version only supports point cloud. 
*/
void Mesh::savePartitionAsOFF(const char *filePath, int edgeId, int *partitionMap) const
{
	vector<Vector3f> coordList;

	for(int i = 0 ; i < this->N ; i++)
	{
		if(partitionMap[i] == edgeId)
			coordList.push_back(vertices[i]->coords);
	}

	FILE *pFile = fopen(filePath, "w");

	fprintf(pFile, "OFF\n");
	fprintf(pFile, "%lu 0 0\n", coordList.size());

	for(int i = 0 ; i < coordList.size() ; i++)
		fprintf(pFile, "%f %f %f\n", coordList[i](0), coordList[i](1), coordList[i](2));

	fclose(pFile);
}

/*
	Parses OFF file (meshFilePath) to construct the mesh. 
*/
void Mesh::parseOFFFile(const char *meshFilePath)
{
	char dummyStr[4];
	int numVertices;
	int numTriangles;
	int numEdges;	// Not used
	int dummyInt;
	FILE *pFile = fopen(meshFilePath, "r");

	fscanf(pFile, "%s", dummyStr);
	fscanf(pFile, "%d %d %d", &numVertices, &numTriangles, &numEdges);

	// Read the vertices.
	for(int i = 0 ; i < numVertices ; i++)
	{
		Vector3f coords;

		fscanf(pFile, "%f %f %f", &coords(0), &coords(1), &coords(2));

		vertices.push_back(new Vertex(i, coords));
	}

	// Read the triangles.
	for(int i = 0 ; i < numTriangles ; i++)
	{
		int v1;
		int v2;
		int v3;
		int edgeid1 = 3 * i;
		int edgeid2 = 3 * i + 1;
		int edgeid3 = 3 * i + 2;

		fscanf(pFile, "%d", &dummyInt);
		fscanf(pFile, "%d %d %d", &v1, &v2, &v3);

		tris.push_back(new Triangle(i, v1, v2, v3));

		edges.push_back(new Edge(edgeid1, v1, v2, -1, i));
		edges.push_back(new Edge(edgeid2, v2, v3, -1, i));
		edges.push_back(new Edge(edgeid3, v3, v1, -1, i));

		vertices[v1]->outgoingEdges.push_back(edgeid1);
		vertices[v2]->outgoingEdges.push_back(edgeid2);
		vertices[v3]->outgoingEdges.push_back(edgeid3);

		avgEdgeLength += (vertices[v1]->coords - vertices[v2]->coords).norm();
		avgEdgeLength += (vertices[v2]->coords - vertices[v3]->coords).norm();
		avgEdgeLength += (vertices[v3]->coords - vertices[v1]->coords).norm();
	}

	fclose(pFile);

	this->N = vertices.size();
	this->M = edges.size();
	avgEdgeLength /= this->M;

	printf("Avg. edge length of the source mesh: %f\n", avgEdgeLength);

	// Set the nornals of the triangles. 
	for(int i = 0 ; i < numTriangles ; i++)
		tris[i]->setNormal(vertices[tris[i]->v1]->coords, vertices[tris[i]->v2]->coords, vertices[tris[i]->v3]->coords);

	for(int i = 0 ; i < this->M ; i++)
	{
		int u = edges[i]->u;
		int v = edges[i]->v;

		for(int j = 0 ; j < vertices[v]->outgoingEdges.size() ; j++)
			if(edges[vertices[v]->outgoingEdges[j]]->v == u)
			{
				edges[i]->opp = edges[vertices[v]->outgoingEdges[j]]->id;

				break;
			}
	}

	// For ray-mesh intersection. 
	for(int i = 0 ; i < tris.size() ; i++)
	{
		tris[i]->vertices = &this->vertices;

		Vector3f P1 = vertices[tris[i]->v1]->coords;
		Vector3f P2 = vertices[tris[i]->v2]->coords;
		Vector3f P3 = vertices[tris[i]->v3]->coords;

		tris[i]->U1 = (P1 - P2).normalized() + (P1 - P3).normalized();
		tris[i]->U2 = (P2 - P3).normalized() + (P2 - P1).normalized();
		tris[i]->U3 = (P3 - P1).normalized() + (P3 - P2).normalized();

		tris[i]->xa_m_xb = P1(0) - P2(0);
		tris[i]->xa_m_xc = P1(0) - P3(0);
		tris[i]->ya_m_yb = P1(1) - P2(1);
		tris[i]->ya_m_yc = P1(1) - P3(1);
		tris[i]->za_m_zb = P1(2) - P2(2);
		tris[i]->za_m_zc = P1(2) - P3(2);
		tris[i]->det_A_term = tris[i]->ya_m_yb * tris[i]->za_m_zc - tris[i]->ya_m_yc * tris[i]->za_m_zb;

		this->triIdxs.push_back(i);

		this->meshArea += ((P1 - P3).cross(P1 - P2)).norm() / 2.0;
	}

	setVerticeNormals();
}

/*
	Corresponds mesh with the IK-skeleton (skel).

	Partitions the mesh into disjoint regions. Every vertex corresponds to a skeleton edge. A mesh 
	vertex corresponds to the closest skeleton edge. The mesh vertices corresponded to the same 
	skeleton edge belongs to the same region. 

	Each element in partitionMap holds id of the edge that element corresponds. e.g., 
	partitionMap[i]=m means i'th mesh vertex corresponds to m'th edge. 

	parametricCoords holds parametric coordinates (For each skeleton edge, u->t=0, v->t=1) of the 
	projection points of the mesh vertices on the IK-skeleton. 
*/
void Mesh::correspondMesh(const Skeleton *skel) const
{
	for(int i = 0 ; i < vertices.size() ; i++)
	{
		tuple<int, float> corrInfo = projectPoint(i, skel);

		int edgeId = get<0>(corrInfo); 
		float paramCoord = get<1>(corrInfo);

		partitions[i] = edgeId;
		parametricCoords[i] = paramCoord;
	}
}

/*
	Helper method for correspondMesh. Projects mesh vertex with id meshVertId onto the IK-skeleton and 
	returns the necessary parameters. Traverses only the edges with even id's. Edges with odd id's are 
	just the opposite half edges. 
*/
tuple<int, float> Mesh::projectPoint(int meshVertId, const Skeleton *skel) const
{
	int edgeId;
	float param;
	float d_min = FLT_MAX;
	Vector3f P = vertices[meshVertId]->coords;

	for(int i = 0 ; i < skel->edges.size() ; i += 2)
	{
		Vector3f projVect;
		Vector3f U = skel->vertices[skel->edges[i]->u]->coords;
		Vector3f V = skel->vertices[skel->edges[i]->v]->coords;
		Vector3f UV = V - U;
		Vector3f UP = P - U;
		float t = UP.dot(UV.normalized()) / UV.norm();

		Vector3f P_ = U + t * UV;	// Projection

		if(t <= 0.0)
			projVect = U - P;
		else if(t >= 1.0)
			projVect = V - P;
		else
			projVect = P_ - P;

		float d_curr = projVect.norm();

		if(d_curr < d_min)
		{
			d_min = d_curr;

			edgeId = i;
			param = t;
		}
	}

	return make_tuple(edgeId, param);
}

/*
	Returns the cartesian coordinates of a point given as parametric coordinates on the curve-skeleton. 
*/
Vector3f Mesh::getCartesianCoords(const Skeleton *skel, int skelEdgeId, float t) const
{
	Vector3f U = skel->vertices[skel->edges[skelEdgeId]->u]->coords;
	Vector3f V = skel->vertices[skel->edges[skelEdgeId]->v]->coords;

	return U + t * (V - U);
}

/*
	Translates point P as much as u (source) of edgeId translated from the source IK to the target IK. 
*/
Vector3f Mesh::applyEdgeTranslation(Vector3f P, int edgeId, const Skeleton *sSkel, const Skeleton *tSkel) const
{
	Vector3f U_s = sSkel->vertices[sSkel->edges[edgeId]->u]->coords;
	Vector3f U_t = tSkel->vertices[tSkel->edges[edgeId]->u]->coords;

	return P + (U_t - U_s);
}

/*
	Rotates P towards the edge so that translated and target edges are aligned. 
*/
Vector3f Mesh::applyEdgeAlignment(Vector3f P, int edgeId, const Skeleton *sSkel, const Skeleton *tSkel) const
{
	Vector3f U_s = sSkel->vertices[sSkel->edges[edgeId]->u]->coords;
	Vector3f U_t = tSkel->vertices[tSkel->edges[edgeId]->u]->coords;
	Vector3f T = U_t - U_s;
	Vector3f V_s_moved = sSkel->vertices[sSkel->edges[edgeId]->v]->coords + T;
	Vector3f V_t = tSkel->vertices[tSkel->edges[edgeId]->v]->coords;
	Vector3f v1 = (V_s_moved - U_t).normalized();
	Vector3f v2 = (V_t - U_t).normalized();
	Vector3f r = (v1.cross(v2)).normalized();

	P -= U_t;

	float theta;
	float v1_dot_v2 = v1.dot(v2);

	if(v1_dot_v2 > 1.0)
		theta = 0.0;
	else if(v1_dot_v2 < -1.0)
		theta = M_PI;
	else
		theta = acos(v1.dot(v2));

	Vector3f Q = rotate(P, theta, r);

	Q += U_t;

	return Q;
}

/*
	Moves point P according to the scaling between the two corresponding edges. 
	TODO: Subject to be removed. 
*/
Vector3f Mesh::applyScaling(Vector3f P, int edgeId, float paramCoord, const Skeleton *sSkel, const Skeleton *tSkel) const
{
	Vector3f U_s = sSkel->vertices[sSkel->edges[edgeId]->u]->coords;
	Vector3f V_s = sSkel->vertices[sSkel->edges[edgeId]->v]->coords;
	Vector3f U_t = tSkel->vertices[tSkel->edges[edgeId]->u]->coords;
	Vector3f V_t = tSkel->vertices[tSkel->edges[edgeId]->v]->coords;
	Vector3f w = (V_t - U_t).normalized();

	if(paramCoord <= 1)
		return P + w * paramCoord * ((V_t - U_t).norm() - (V_s - U_s).norm());
	else
		return P + w * (paramCoord - 1.0) * ((V_t - U_t).norm() - (V_s - U_s).norm());
}

void Mesh::printBorderEdgeLengthDiff(const Mesh *sMesh, int *partitionMap)
{
	bool visitedEdges[edges.size()];
	float totalBorderEdgeInc = 0.0;
	float sourceBorderEdgeLen = 0.0;
	float targetBorderEdgeLen = 0.0;
	vector<int> borderEdges;

	for(int i = 0 ; i < edges.size() ; i++)
		visitedEdges[i] = false;

	for(int i = 0 ; i < edges.size() ; i++)
	{
		if(partitionMap[edges[i]->u] != partitionMap[edges[i]->v] && !visitedEdges[i])
		{
			borderEdges.push_back(i);
			visitedEdges[i] = true;
			visitedEdges[edges[i]->opp] = true;
		}
	}

	for(int i = 0 ; i < borderEdges.size() ; i++)
	{
		Vector3f U_s = sMesh->vertices[sMesh->edges[borderEdges[i]]->u]->coords;
		Vector3f V_s = sMesh->vertices[sMesh->edges[borderEdges[i]]->v]->coords;
		Vector3f U_t = this->vertices[this->edges[borderEdges[i]]->u]->coords;
		Vector3f V_t = this->vertices[this->edges[borderEdges[i]]->v]->coords;

		float l_s = (U_s - V_s).norm();
		float l_t = (U_t - V_t).norm();

		sourceBorderEdgeLen += l_s;
		targetBorderEdgeLen += l_t;
		totalBorderEdgeInc += l_t - l_s;
	}

	printf("%f %f %f\n", sourceBorderEdgeLen, targetBorderEdgeLen, totalBorderEdgeInc);
}

/*
	TODO: Only version 1.1. Will be improved. 
*/
void Mesh::markFreeVertices(float ratio, bool *freeVertices, const Skeleton *skel) const
{
	// Only a subset of the initial mesh vertices are free vertices. 
	for(int i = 0 ; i < this->N ; i++)
	{
		int corrEdge;

		freeVertices[i] = false;	// Initialize as handle. 
		vertices[i]->region = 0;	// 0 means rigid region. 

		if(skel->edges[partitions[i]]->length < skel->avgEdgeLength)
		{
			if(!skel->vertices[skel->edges[partitions[i]]->u]->isTerminal() && !skel->vertices[skel->edges[partitions[i]]->v]->isTerminal())
			{
				freeVertices[i] = true;
				vertices[i]->region = -1;	// -1 means not assigned. 

				continue;
			}
		}

		if(parametricCoords[i] < ratio)
			corrEdge = partitions[i];
		else if(parametricCoords[i] > 1.0 - ratio)
			corrEdge = partitions[i] + 1;
		else
			continue;	// Leave it as handle. 

		if(!skel->vertices[skel->edges[corrEdge]->u]->isTerminal())
		{
			freeVertices[i] = true;
			vertices[i]->region = -1;
		}
	}

	// All the support vertices are handles. 
	for(int i = this->N ; i < vertices.size() ; i++)
	{
		freeVertices[i] = false;
		vertices[i]->region = 0;
	}
}

/*
	Computes weight for the projection edge between meshVertId and meshVertId + N
*/
float Mesh::computeProjEdgeWeight(int meshVertId) const
{
	unsigned long N = vertices.size();
	unsigned long E = edges.size();
	float w_ij = 0.0;
	Vector3f P_i = vertices[meshVertId]->coords;
	Vector3f P_s = vertices[meshVertId + N]->coords;
	float l_is = (P_i - P_s).norm();

	for(int j = 0 ; j < vertices[meshVertId]->outgoingEdges.size() ; j++)
	{
		if(vertices[meshVertId]->outgoingEdges[j] < E)
		{
			Vector3f P_j = vertices[edges[vertices[meshVertId]->outgoingEdges[j]]->v]->coords;
			float l_ij = (P_i - P_j).norm();

			w_ij += l_ij / l_is;
		}
	}

	return w_ij;
}

void Mesh::saveNonRigidPartsAsOFF(const char *filePath, bool *freeVertices) const
{
	vector<Vector3f> coordList;

	for(int i = 0 ; i < vertices.size() ; i++)
	{
		if(!freeVertices[i])
			coordList.push_back(vertices[i]->coords);
	}

	FILE *pFile = fopen(filePath, "w");

	fprintf(pFile, "OFF\n");
	fprintf(pFile, "%lu 0 0\n", coordList.size());

	for(int i = 0 ; i < coordList.size() ; i++)
		fprintf(pFile, "%f %f %f\n", coordList[i](0), coordList[i](1), coordList[i](2));

	fclose(pFile);
}

/*
	Computes cot(alpha) where alpha is the opposite angle of edge with id edgeId. 
	TODO: Speed up by holding edge lengths from the beginning. 
*/
float Mesh::computeCotan(int edgeId) const
{
	int vert1 = edges[edgeId]->u;
	int vert2 = edges[edgeId]->v;
	int vert3 = tris[edges[edgeId]->faceId]->getThirdVertex(vert1, vert2);

	float a = (vertices[vert1]->coords - vertices[vert2]->coords).norm();
	float b = (vertices[vert2]->coords - vertices[vert3]->coords).norm();
	float c = (vertices[vert3]->coords - vertices[vert1]->coords).norm();

	float cos_alpha = (b * b + c * c - a * a) / (2 * b * c);
	float sin_alpha = sqrt(1 - cos_alpha * cos_alpha);

	return cos_alpha / sin_alpha;
}

/*
	Computes optimum rotation matrix mapping patch i of the source mesh to the patch i of the target mesh. 
	P is the previous point set. 
*/
void Mesh::computeOptimumRotation(const Mesh *sMesh, int idx, SparseMatrix<float> & W, vector<vector<int>> & neighbourhoodInfo, Matrix3f & R) const
{
	int n = neighbourhoodInfo[idx].size();	// Number of neighbours. 
	MatrixXf X(3, n);	// 3d Source patch. 
	MatrixXf Y(3, n);	// 3d Target patch. 
	MatrixXf D = MatrixXf::Zero(n, n);	// Diagonal matrix containing w_ij's. 

	Vector3f p_i = sMesh->vertices[idx]->coords;
	Vector3f p_i_ = this->vertices[idx]->coords;

	// Arrange X, Y and D. 
	for(int i = 0 ; i < n ; i++)
	{
		Vector3f p_j = sMesh->vertices[neighbourhoodInfo[idx][i]]->coords;
		Vector3f p_j_ = this->vertices[neighbourhoodInfo[idx][i]]->coords;

		X.col(i) = p_i - p_j;
		Y.col(i) = p_i_ - p_j_;
		D(i, i) = W.coeffRef(idx, neighbourhoodInfo[idx][i]);
	}

	const Matrix<float, 3, 3> S = X * D * Y.transpose();
	Matrix<float, 3, 3> U;
	Matrix<float, 3, 1> Sigma;
	Matrix<float, 3, 3> V;

	wunderSVD3x3(S, U, Sigma, V);

	R = V * U.transpose();
}

/*
	Computes optimum smooth local rotations (SR-ARAP).  
*/
Matrix3f Mesh::computeOptimumSmoothRotation(const Mesh *sMesh, int idx, SparseMatrix<float> & W, vector<Matrix3f> & Rs) const
{
	float alpha = 0.01;
	Vector3f P_i = sMesh->vertices[idx]->coords;
	Vector3f Q_i = this->vertices[idx]->coords;
	MatrixXf S = MatrixXf::Zero(3, 3);
	MatrixXf R(3, 3);
	JacobiSVD<MatrixXf> svd;

	for(int j = 0 ; j < vertices[idx]->outgoingEdges.size() ; j++)
	{
		int neighIdx = edges[vertices[idx]->outgoingEdges[j]]->v;
		Vector3f P_j = sMesh->vertices[neighIdx]->coords;
		Vector3f Q_j = this->vertices[neighIdx]->coords;

		S += W.coeffRef(idx, neighIdx) * (P_j - P_i) * (Q_j - Q_i).transpose();
	}

	S *= 2;

	for(int j = 0 ; j < this->vertices[idx]->outgoingEdges.size() ; j++)
	{
		int neighIdx = this->edges[vertices[idx]->outgoingEdges[j]]->v;

		S += 4 * alpha * sMesh->meshArea * Rs[neighIdx].transpose();
	}
	
	svd.compute(S, Eigen::ComputeThinU | Eigen::ComputeThinV);
	if(!svd.computeU() || !svd.computeV())
	{
		printf("Decomposition error.\n");
		exit(0);
	}

	MatrixXf V = svd.matrixV();
	MatrixXf U = svd.matrixU();
	MatrixXf VUt = V * U.transpose();
	MatrixXf M(3, 3);

	M.setIdentity();

	M(2, 2) = VUt.determinant();

	R = V * M * U.transpose();

	return R;
}

/*
	Updates vertex positions in the previous step. 
*/
void Mesh::updatePositions(const Mesh *sMesh, vector<Matrix3f> & R, SparseLU<SparseMatrix<float>, COLAMDOrdering<int>> & solver, SparseMatrix<float> & W, SparseMatrix<float> & L, bool *freeVertices) const
{
	VectorXf b_x(sMesh->vertices.size());
	VectorXf b_y(sMesh->vertices.size());
	VectorXf b_z(sMesh->vertices.size());
	VectorXf P_x_(sMesh->vertices.size());
	VectorXf P_y_(sMesh->vertices.size());
	VectorXf P_z_(sMesh->vertices.size());

	// Computation of b_x, b_y, and b_z. 
	for(int i = 0 ; i < sMesh->vertices.size() ; i++)
	{
		if(!freeVertices[i])	// i is handle. 
		{
			b_x(i) = vertices[i]->coords(0);	// x coordinates. 
			b_y(i) = vertices[i]->coords(1);	// y coordiantes. 
			b_z(i) = vertices[i]->coords(2);	// z coordinates. 

			continue;
		}

		// i is not handle. 
		Vector3f b_i(0.0, 0.0, 0.0);
		Vector3f p_i = sMesh->vertices[i]->coords;

		for(int j = 0 ; j < sMesh->vertices[i]->outgoingEdges.size() ; j++)
		{
			int neighInd = sMesh->edges[sMesh->vertices[i]->outgoingEdges[j]]->v;
			Vector3f p_j = sMesh->vertices[neighInd]->coords;

			b_i += 0.5 * W.coeffRef(i, neighInd) * ((R[i] + R[neighInd]) * (p_i - p_j));
		}

		b_x(i) = b_i(0);	// x coordinates. 
		b_y(i) = b_i(1);	// y coordinates. 
		b_z(i) = b_i(2);	// z coordinates. 
	}

	thread th1(&Mesh::solveForPositions, this, std::ref(b_x), std::ref(solver), std::ref(P_x_));	// Solve for L*p_x=b_x.
	thread th2(&Mesh::solveForPositions, this, std::ref(b_y), std::ref(solver), std::ref(P_y_));	// Solve for L*p_y=b_y.
	thread th3(&Mesh::solveForPositions, this, std::ref(b_z), std::ref(solver), std::ref(P_z_));	// Solve for L*p_z=b_z. 

	th1.join();
	th2.join();
	th3.join();

	// Arrange P_.
	for(int i = 0 ; i < sMesh->vertices.size() ; i++)
	{
		this->vertices[i]->coords(0) = P_x_[i];
		this->vertices[i]->coords(1) = P_y_[i];
		this->vertices[i]->coords(2) = P_z_[i];
	}
}

/*
	Solves for P in system Ax=b. 
*/
void Mesh::solveForPositions(VectorXf & b, SparseLU<SparseMatrix<float>, COLAMDOrdering<int>> & solver, VectorXf & x) const
{
	x = solver.solve(b);
}

float Mesh::computeError(const Mesh *sMesh, vector<Matrix3f> & R, SparseMatrix<float> & W, bool onlyMesh) const
{
	float totalError = 0.0;
	int border = onlyMesh ? this->N : vertices.size();

	for(int i = 0 ; i < border ; i++)
		totalError += computePerCellError(sMesh, i, R[i], W, onlyMesh);

	return totalError;
}

float Mesh::computePerCellError(const Mesh *sMesh, int vert, Matrix3f & R, SparseMatrix<float> & W, bool onlyMesh) const
{
	float perCellError = 0.0;

	for(int j = 0 ; j < vertices[vert]->outgoingEdges.size() ; j++)
	{
		int neighIdx = edges[vertices[vert]->outgoingEdges[j]]->v;
		Vector3f P_j = sMesh->vertices[neighIdx]->coords;
		Vector3f P_j_ = vertices[neighIdx]->coords;

		if(!onlyMesh || neighIdx < this->N)
		{
			// TODO: What if w_ij is negative? |w_ij| is used here. 
			perCellError += fabs(W.coeffRef(vert, neighIdx)) * pow(((this->vertices[vert]->coords - P_j_) - R * (sMesh->vertices[vert]->coords - P_j)).norm(), 2);
		}
	}

	return perCellError;
}

void Mesh::taubinMeshFairingOneStep(void)
{
	float lambda = 0.6307;
	float mu = -0.6732;
	Vector3f regulators[this->N];

	for(int i = 0 ; i < this->N ; i++)	// Process only mesh vertices. 
	{
		float dist;
		float totalWeight = 0.0;
		float weight;
		Vector3f disp;
		Vector3f P_i = this->vertices[i]->coords;

		regulators[i] = Vector3f::Zero();

		for(int j = 0 ; j < this->vertices[i]->outgoingEdges.size() ; j++)
		{
			int edgeIdx = this->vertices[i]->outgoingEdges[j];

			if(edgeIdx < this->M)
			{
				Vector3f P_j = this->vertices[this->edges[edgeIdx]->v]->coords;

				disp = P_j - P_i;
				dist = disp.norm();

				if(dist == 0.0)
				{
					regulators[i] = Vector3f::Zero();

					break;
				}

				weight = 1.0 / dist;
				totalWeight += weight;
				disp *= weight;
				regulators[i] += disp;
			}
		}

		regulators[i] /= totalWeight;
	}

	for(int i = 0 ; i < this->N ; i++)
		this->vertices[i]->coords += lambda * regulators[i];

	for(int i = 0 ; i < this->N ; i++)
	{
		float dist;
		float totalWeight = 0.0;
		float weight;
		Vector3f disp;
		Vector3f P_i = this->vertices[i]->coords;

		regulators[i] = Vector3f::Zero();

		for(int j = 0 ; j < this->vertices[i]->outgoingEdges.size() ; j++)
		{
			int edgeIdx = this->vertices[i]->outgoingEdges[j];

			if(edgeIdx < this->M)
			{
				Vector3f P_j = this->vertices[this->edges[edgeIdx]->v]->coords;

				disp = P_j - P_i;
				dist = disp.norm();

				if(dist == 0.0)
				{
					regulators[i] = Vector3f::Zero();

					break;
				}

				weight = 1.0 / dist;
				totalWeight += weight;
				disp *= weight;
				regulators[i] += disp;
			} 
		}

		regulators[i] /= totalWeight;
	}

	for(int i = 0 ; i < this->N ; i++)
		this->vertices[i]->coords += mu * regulators[i];
}

/*
	Computes the initial solution for the transferred mesh. 
*/
void Mesh::computeInitialSolution(Mesh *sMesh, const Skeleton *sSkel, const Skeleton *tSkel)
{
	vector<vector<float>> w_LBS;

	//w_LBS = sMesh->computeHeatBoneWeights(sSkel);
	loadCoeffs("../tr_reg_039_stillIm.coeffs", w_LBS);

	vector<Matrix4f> T;	// 4x4 matrices associated with bones (translation+rotation). 

	for(int i = 0 ; i < sSkel->edges.size() ; i += 2)
	{
		Matrix4f T_i;

		//if(i == 2)
		//	T_i = computeTransformationMatrix(sSkel, tSkel, M_PI, i);
		//else
			T_i = computeTransformationMatrix(sSkel, tSkel, 0.0, i);

		T.push_back(T_i);
		T.push_back(T_i);
	}

	for(int i = 0 ; i < this->N ; i++)
	{
		Vector4f sPoint(sMesh->vertices[i]->coords(0), sMesh->vertices[i]->coords(1), sMesh->vertices[i]->coords(2), 1.0);
		Vector4f tPoint(0.0, 0.0, 0.0, 0.0);

		for(int j = 0 ; j < w_LBS[i].size() ; j += 2)
			tPoint += w_LBS[i][j] * T[j] * sPoint;

		this->vertices[i]->coords(0) = tPoint(0);
		this->vertices[i]->coords(1) = tPoint(1);
		this->vertices[i]->coords(2) = tPoint(2);
	}
}

/*
	Computes the transformation matrix mapping the bone with edgeId in sSkel to tSkel. 
*/
Matrix4f Mesh::computeTransformationMatrix(const Skeleton *sSkel, const Skeleton *tSkel, float phi, int edgeId)
{
	Matrix4f T = MatrixXf::Identity(4, 4);
	
	// Compute translation. 
	Vector3f U_s = sSkel->vertices[sSkel->edges[edgeId]->u]->coords;
	Vector3f U_t = tSkel->vertices[tSkel->edges[edgeId]->u]->coords;

	// Compute rotation. 
	Vector3f V_s_moved = sSkel->vertices[sSkel->edges[edgeId]->v]->coords + (U_t - U_s);
	Vector3f V_t = tSkel->vertices[tSkel->edges[edgeId]->v]->coords;
	Vector3f v1 = (V_s_moved - U_t).normalized();
	Vector3f v2 = (V_t - U_t).normalized();
	Vector3f r = (v1.cross(v2)).normalized();

	float costheta = v1.dot(v2);

	if(costheta > 1.0)
		costheta = 1.0;
	else if(costheta < -1.0)
		costheta = -1.0;

	float sintheta = sqrt(1.0 - costheta * costheta);

	T(0, 0) = costheta + (1 - costheta) * r(0) * r(0);
	T(0, 1) = (1 - costheta) * r(0) * r(1) - r(2) * sintheta;
	T(0, 2) = (1 - costheta) * r(0) * r(2) + r(1) * sintheta;

	T(1, 0) = (1 - costheta) * r(0) * r(1) + r(2) * sintheta;
	T(1, 1) = costheta + (1 - costheta) * r(1) * r(1);
	T(1, 2) = (1 - costheta) * r(1) * r(2) - r(0) * sintheta;

	T(2, 0) = (1 - costheta) * r(0) * r(2) - r(1) * sintheta;
	T(2, 1) = (1 - costheta) * r(1) * r(2) + r(0) * sintheta;
	T(2, 2) = costheta + (1 - costheta) * r(2) * r(2);

	// Rotation around edge.
	Vector3f V_s = sSkel->vertices[sSkel->edges[edgeId]->v]->coords;
	Vector3f r_e = (V_s - U_s).normalized();

	float cosphi = cos(phi);

	if(cosphi > 1.0)
		cosphi = 1.0;
	else if(cosphi < -1.0)
		cosphi = -1.0;

	float sinphi = sqrt(1.0 - cosphi * cosphi);

	Matrix4f T_e = MatrixXf::Identity(4, 4);

	T_e(0, 0) = cosphi + (1 - cosphi) * r_e(0) * r_e(0);
   	T_e(0, 1) = (1 - cosphi) * r_e(0) * r_e(1) - r_e(2) * sinphi;
	T_e(0, 2) = (1 - cosphi) * r_e(0) * r_e(2) + r_e(1) * sinphi;

	T_e(1, 0) = (1 - cosphi) * r_e(0) * r_e(1) + r_e(2) * sinphi;
	T_e(1, 1) = cosphi + (1 - cosphi) * r_e(1) * r_e(1);
	T_e(1, 2) = (1 - cosphi) * r_e(1) * r_e(2) - r_e(0) * sinphi;

	T_e(2, 0) = (1 - cosphi) * r_e(0) * r_e(2) - r_e(1) * sinphi;
	T_e(2, 1) = (1 - cosphi) * r_e(1) * r_e(2) + r_e(0) * sinphi;
	T_e(2, 2) = cosphi + (1 - cosphi) * r_e(2) * r_e(2);

	Matrix4f T_pre = MatrixXf::Identity(4, 4);

	T_pre(0, 3) = -U_s(0);
	T_pre(1, 3) = -U_s(1);
	T_pre(2, 3) = -U_s(2);

	Matrix4f T_pst = MatrixXf::Identity(4, 4);

	T_pst(0, 3) = U_t(0);
	T_pst(1, 3) = U_t(1);
	T_pst(2, 3) = U_t(2);

	return T_pst * T * T_e * T_pre;
}

void Mesh::loadCoeffs(const char *filePath, vector<vector<float>> & coeffs)
{
	int numVerts;
	int numBones;
	FILE *pFile = fopen(filePath, "r");

	fscanf(pFile, "%d%d", &numVerts, &numBones);

	coeffs.resize(numVerts, vector<float>(numBones));

	for(int i = 0 ; i < numVerts ; i++)
		for(int j = 0 ; j < numBones ; j++)
			fscanf(pFile, "%f", &coeffs[i][j]);

	fclose(pFile);
}
