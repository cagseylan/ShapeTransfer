#include "Skeleton.h"

Skeleton::Skeleton()
{}

/*
	Skeleton construction by reading .cg file. 
*/
Skeleton::Skeleton(const char *filePath)
{
	char currChar;
	int currVertexCnt = 0;
	int currEdgeCnt = 0;
	FILE *pFile = fopen(filePath, "r");

	// Parse the .cg file. 
	while(fscanf(pFile, "%c", &currChar) != EOF)
	{
		if(currChar == EOF)
			break;
		switch(currChar)
		{
			case '#':	// Only comment, do nothing.
				while(1)
				{
					fscanf(pFile, "%c", &currChar);
					if(currChar == '\n')
						break;
				}
				break;
			case 'v':	// Vertex.
			{
				Vector3f coords;

				fscanf(pFile, "%f %f %f", &coords(0), &coords(1), &coords(2));
				vertices.push_back(new Vertex(currVertexCnt, coords));
				currVertexCnt++;

				break;
			}
			case 'e':	// Edge.
			{
				int v1;
				int v2;

				fscanf(pFile, "%d %d", &v1, &v2);
				edges.push_back(new Edge(2 * currEdgeCnt, v1 - 1, v2 - 1, 2 * currEdgeCnt + 1, -1));
				edges.push_back(new Edge(2 * currEdgeCnt + 1, v2 - 1, v1 - 1, 2 * currEdgeCnt, -1));
				currEdgeCnt++;

				break;
			}
			default:
				break;
		}
	}

	for(int i = 0 ; i < edges.size() ; i++)
		vertices[edges[i]->u]->outgoingEdges.push_back(edges[i]->id);

	fclose(pFile);

	for(int i = 0 ; i < edges.size() ; i += 2)
	{
		Vector3f U = vertices[edges[i]->u]->coords;
		Vector3f V = vertices[edges[i]->v]->coords;

		float edgeLength = (V - U).norm();

		edges[i]->length = edgeLength;
		edges[i + 1]->length = edgeLength;
		avgEdgeLength += edgeLength;
	}

	avgEdgeLength /= (edges.size() / 2.0);
}

/*
	Saves the skeleton in CG format into the path fileName.
*/
void Skeleton::saveAsCG(const char *fileName)
{
	FILE *pFile = fopen(fileName, "w");

	fprintf(pFile, "# D:3 ");
	fprintf(pFile, "NV:%lu ", vertices.size());
	fprintf(pFile, "NE:%lu\n", edges.size() / 2);

	for(int i = 0 ; i < vertices.size() ; i++)
		fprintf(pFile, "v %f %f %f\n", vertices[i]->coords(0), vertices[i]->coords(1), vertices[i]->coords(2));
	for(int i = 0 ; i < edges.size() ; i += 2)
		fprintf(pFile, "e %d %d\n", edges[i]->u + 1, edges[i]->v + 1);

	fclose(pFile);
}

/*
	Find index of trunk (The edge with largest number of outgoing edges). 
	If there is more than one, just get the first one. 
*/
int Skeleton::findTrunk(void) const
{
	int trunkIdx;
	int degTrunk = 0;

	for(int i = 0 ; i < edges.size() ; i += 2)
	{
		int u = edges[i]->u;
		int v = edges[i]->v;
		int totDeg = 0;

		int deg_u = vertices[u]->outgoingEdges.size();
		int deg_v = vertices[v]->outgoingEdges.size();

		totDeg = deg_u + deg_v - 2;

		if(totDeg > degTrunk)
		{
			degTrunk = totDeg;
			trunkIdx = i;
		}
	}

	return trunkIdx;
}

void Skeleton::getSkeletonVertices(MatrixXd & skeletonVertices) const
{
	skeletonVertices.resize(vertices.size(), 3);

	for(int i = 0 ; i < vertices.size() ; i++)
	{
		skeletonVertices(i, 0) = vertices[i]->coords(0);
		skeletonVertices(i, 1) = vertices[i]->coords(1);
		skeletonVertices(i, 2) = vertices[i]->coords(2);
	}
}

vector<tuple<int, int>> Skeleton::getSkeletonEdges(void) const
{
	vector<tuple<int, int>> skeletonEdges;

	for(int i = 0 ; i < edges.size() ; i += 2)
		skeletonEdges.push_back(make_tuple(edges[i]->u, edges[i]->v));

	return skeletonEdges;
}

/*
	Samples the curve-skeleton with a gap of stepSize. 
	Returns vector of tuples where the first element of each tuple is the coordinates of the sample point, 
	and the second element is the edge direction at the sample point. 
*/
vector<tuple<Vector3f, Vector3f>> Skeleton::sampleSkeleton(float stepSize) const
{
	vector<tuple<Vector3f, Vector3f>> skeletonSamples;
	vector<vector<int>> segments = computeSegments();

	for(int i = 0 ; i < segments.size() ; i++)
	{
		vector<tuple<Vector3f, Vector3f>> segmentSamples = sampleSegment(segments[i], stepSize);

		skeletonSamples.insert(skeletonSamples.end(), segmentSamples.begin(), segmentSamples.end());
	}

	return skeletonSamples;
}

/*
	stepSize: gap between samples in the curve-skeleton. 
	rayPerSample: Number of rays emitted from each sample. 
	Returns vector of tuples, where first element in each tuple is the ray origin, and 
	the vector in each tuple holds ray directions. 
*/
vector<tuple<Vector3f, vector<Vector3f>>> Skeleton::computeRays(float stepSize, int raysPerSample) const
{
	vector<tuple<Vector3f, vector<Vector3f>>> rays;
	vector<tuple<Vector3f, Vector3f>> samples = sampleSkeleton(stepSize);

	for(int i = 0 ; i < samples.size() ; i++)
	{
		Vector3f org = get<0>(samples[i]);	// Origin of the i'th ray. 
		Vector3f edgeDirection = get<1>(samples[i]); // Direction of the edge at org. 

		vector<Vector3f> rayDirections = computeDirections(org, edgeDirection, raysPerSample);	// Directions of the i'th ray. 

		rays.push_back(make_tuple(org, rayDirections));
	}

	return rays;
}

/*
	Returns length of the edge with edgeId. 
*/
float Skeleton::edgeLength(int edgeId)
{
	Vector3f u = vertices[edges[edgeId]->u]->coords;
	Vector3f v = vertices[edges[edgeId]->v]->coords;

	return (v - u).norm();
}

/*
	Computes the segments of the skeleton in a list of list. Inner list contains vertex ids of vertices in the segment in order. 
*/
vector<vector<int>> Skeleton::computeSegments(void) const
{
	bool visitedEdges[edges.size()];
	vector<vector<int>> segments;

	for(int i = 0 ; i < edges.size() ; i++)
		visitedEdges[i] = false;

	for(int i = 0 ; i < vertices.size() ; i++)
	{
		if(!vertices[i]->isRegular())	// Terminal or joint vertex. 
		{
			for(int j = 0 ; j < (vertices[i]->outgoingEdges).size() ; j++)
			{
				int currEdge;
				int currVertex;
				vector<int> skeletonSegment;

				currEdge = vertices[i]->outgoingEdges[j];
				
				if(visitedEdges[currEdge])
					continue;

				currVertex = vertices[i]->id;
				skeletonSegment.push_back(currVertex);

				currVertex = edges[currEdge]->v;

				visitedEdges[currEdge] = true;
				visitedEdges[edges[currEdge]->opp] = true;

				while(vertices[currVertex]->isRegular())
				{
					skeletonSegment.push_back(currVertex);

					if(vertices[currVertex]->outgoingEdges[0] == edges[currEdge]->opp)
						currEdge = vertices[currVertex]->outgoingEdges[1];
					else
						currEdge = vertices[currVertex]->outgoingEdges[0];

					currVertex = edges[currEdge]->v;

					visitedEdges[currEdge] = true;
					visitedEdges[edges[currEdge]->opp] = true;
				}

				skeletonSegment.push_back(currVertex);
				segments.push_back(skeletonSegment);
			}
		}
	}

	return segments;
}

/*
	Computes distance between a line segment (UV) and a point (P). 
*/
double Skeleton::computeDistance2LineSeg(const Vector3f & U, const Vector3f & V, const Vector3f & P) const
{
	Vector3f UV = V - U;
	Vector3f UP = P - U;
	double t = UV.dot(UP) / UV.dot(UV);

	if(t <= 0.0)
		return (P - U).norm();
	else if(t >= 1.0)
		return (P - V).norm();
	else
		return (P - (U + t * UV)).norm();
}

/*
	Samples a skeleton segment (which is an arc) with a spacing of stepSize.
	Returns a vector of tuples where the first element of each tuple is the coordinates of the sample, 
	and the second element is the edge direction at the sample point. 
*/
vector<tuple<Vector3f, Vector3f>> Skeleton::sampleSegment(vector<int> segment, float stepSize) const
{

	Vector3f p_t;
	vector<tuple<Vector3f, Vector3f>> samples;

	p_t = vertices[segment[0]]->coords;

	int i = 1;	// Index of the next vertex on the segment. 
	if(vertices[segment[0]]->isTerminal())
	{
		Vector3f edgeDir = (vertices[segment[i]]->coords - p_t).normalized();

		samples.push_back(make_tuple(p_t, edgeDir));
	}

	while(true)
	{
		// Move p_t on the arc as step size.
		float dist2Next = (vertices[segment[i]]->coords - p_t).norm();

		if(stepSize < dist2Next)
		{
			Vector3f dir = (vertices[segment[i]]->coords - p_t).normalized();

			p_t = p_t + dir * stepSize;
			samples.push_back(make_tuple(p_t, dir));
		}
		else if(stepSize == dist2Next)
		{
			Vector3f dir(0, 0, 0);

			p_t = vertices[segment[i]]->coords;
			dir += (p_t - vertices[segment[i - 1]]->coords).normalized();
			i++;

			if(i == segment.size())
			{
				samples.push_back(make_tuple(p_t, dir));

				break;
			}

			dir += (vertices[segment[i]]->coords - p_t).normalized();
			samples.push_back(make_tuple(p_t, dir.normalized()));
		}
		else	// stepSize >= dist2Next
		{
			Vector3f dir;

			p_t = vertices[segment[i]]->coords;
			i++;

			if(i == segment.size())
			{
				dir = (p_t - vertices[segment[i - 2]]->coords).normalized();
				samples.push_back(make_tuple(p_t, dir));

				break;
			}

			double residual = stepSize - dist2Next;
			dir = (vertices[segment[i]]->coords - p_t).normalized();

			p_t = p_t + dir * residual;
			samples.push_back(make_tuple(p_t, dir));
		}
	}

	return samples;
}

/*
	Computes directions of the rays emanating from origin with edgeDirection which is the direction of the edge at origin. 
*/
vector<Vector3f> Skeleton::computeDirections(Vector3f origin, Vector3f edgeDirection, int raysPerSample) const
{
	float angleStep = 2 * M_PI / raysPerSample;
	vector<Vector3f> directions;

	Vector3f P(10.0, 10.0, 10.0);	// Dummy point which will help us compute the seed ray perpendicular to edgeDirection.
	Vector3f P_ = origin + edgeDirection * (P - origin).dot(edgeDirection);	// Projection of P onto the ray. 
	Vector3f n = (P - P_).normalized();	// Seed orthogonal vector. 
	Vector3f r = edgeDirection;

	directions.push_back(n);

	for(int i = 1 ; i < raysPerSample ; i++)
	{
		float costheta = cos(angleStep * i);
		float sintheta = sin(angleStep * i);
		Vector3f n_(0.0, 0.0, 0.0);

		n_(0) += (costheta + (1 - costheta) * r(0) * r(0)) * n(0);
   		n_(0) += ((1 - costheta) * r(0) * r(1) - r(2) * sintheta) * n(1);
		n_(0) += ((1 - costheta) * r(0) * r(2) + r(1) * sintheta) * n(2);

		n_(1) += ((1 - costheta) * r(0) * r(1) + r(2) * sintheta) * n(0);
		n_(1) += (costheta + (1 - costheta) * r(1) * r(1)) * n(1);
		n_(1) += ((1 - costheta) * r(1) * r(2) - r(0) * sintheta) * n(2);

		n_(2) += ((1 - costheta) * r(0) * r(2) - r(1) * sintheta) * n(0);
		n_(2) += ((1 - costheta) * r(1) * r(2) + r(0) * sintheta) * n(1);
		n_(2) += (costheta + (1 - costheta) * r(2) * r(2)) * n(2);

		n_.normalize();

		directions.push_back(n_);
	}

	return directions;
}
