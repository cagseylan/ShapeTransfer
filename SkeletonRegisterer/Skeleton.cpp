#include "Skeleton.h"

Skeleton::Skeleton(void)
{}

/*
	Constructs skeleton by reading .cg file.
*/
Skeleton::Skeleton(const char *fileName)
{
	char currChar;
	int currVertexCnt = 0;
	int currEdgeCnt = 0;
	double minEdgeLength = FLT_MAX;
	FILE *pFile = fopen(fileName, "r");

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
				double edgeLength;

				fscanf(pFile, "%d %d", &v1, &v2);
				edgeLength = (vertices[v1 - 1]->coords - vertices[v2 - 1]->coords).norm();
				totalLength += edgeLength;
				edges.push_back(new Edge(2 * currEdgeCnt, v1 - 1, v2 - 1, 2 * currEdgeCnt + 1, edgeLength));
				edges.push_back(new Edge(2 * currEdgeCnt + 1, v2 - 1, v1 - 1, 2 * currEdgeCnt, edgeLength));
				currEdgeCnt++;

				if(edgeLength < minEdgeLength)
				{
					minEdgeLength = edgeLength;
					shortestEdgeIdx = 2 * currEdgeCnt;
				}

				break;
			}
			default:
				break;
		}
	}

	for(int i = 0 ; i < edges.size() ; i++)
		vertices[edges[i]->u]->outgoingEdges.push_back(edges[i]->id);

	fclose(pFile);

	computeSegments();
	computeBranchInfo();
}

Skeleton::Skeleton(const Skeleton *embSkel, const Skeleton *curveSkel, vector<int> pairs)
{
	this->branchInfo = embSkel->branchInfo;
	this->segments = embSkel->segments;

	for(int i = 0 ; i < embSkel->vertices.size() ; i++)
	{
		this->vertices.push_back(new Vertex(i, curveSkel->vertices[pairs[i]]->coords));
		this->vertices[i]->outgoingEdges = embSkel->vertices[i]->outgoingEdges;
	}

	for(int i = 0 ; i < embSkel->edges.size() ; i++)
	{
		this->edges.push_back(new Edge(embSkel->edges[i]->id, embSkel->edges[i]->u, embSkel->edges[i]->v, embSkel->edges[i]->opp, 0.0));
		this->edges[i]->length = (this->vertices[this->edges[i]->u]->coords - this->vertices[this->edges[i]->v]->coords).norm();
		this->totalLength += this->edges[i]->length;
	}
}

/*
	Saves the skeleton in OFF format. 
*/
void Skeleton::saveAsOFF(const char *fileName)
{
	FILE *pFile = fopen(fileName, "w");

	fprintf(pFile, "OFF\n");
	fprintf(pFile, "%lu 0 0\n", vertices.size());

	for(int i = 0 ; i < vertices.size() ; i++)
		fprintf(pFile, "%f %f %f\n", vertices[i]->coords(0), vertices[i]->coords(1), vertices[i]->coords(2));

	fclose(pFile);
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
	Computes the dissimilarity between two curve-skeletons. The lower the better. 
	Basically, it does integration over curve skeletons. 
*/
double Skeleton::computeDissimilarity(const Skeleton *skel) const
{
	double minEdgeLength1 = this->edges[this->shortestEdgeIdx]->length;
	double minEdgeLength2 = skel->edges[skel->shortestEdgeIdx]->length;
	double stepSize = (minEdgeLength1 < minEdgeLength2 ? minEdgeLength1 : minEdgeLength2) / 4.0;

	return this->computeAvgDistance(stepSize, skel) + skel->computeAvgDistance(stepSize, this);
}

/*
	Computes Hausdorff distance between two curve-skeletons. The lower the better.
*/
double Skeleton::computeHausdorffDistance(const Skeleton *skel) const
{
	double minEdgeLength1 = this->edges[this->shortestEdgeIdx]->length;
	double minEdgeLength2 = skel->edges[skel->shortestEdgeIdx]->length;
	double stepSize = (minEdgeLength1 < minEdgeLength2 ? minEdgeLength1 : minEdgeLength2) / 4.0;

	double firstComponent = this->supInt(stepSize, skel);
	double secondComponent = skel->supInt(stepSize, this);

	return (firstComponent > secondComponent ? firstComponent : secondComponent);
}

/*
	Computes shortest path between source and dest. Returns indexes of the vertices on the path in a sorted manner. 
	Writes the shortest distances of each vertex into distanceInfo.
	TODO: May need to handle cycles for shapes with more than one genus. Then, we need Dijkstra's shortest path algorithm. 
*/
vector<int> Skeleton::computeShortestPath_BFS(int source, int dest, double *distanceInfo) const
{
	bool visited[vertices.size()];
	int parentInfo[vertices.size()];
	vector<int> path;
	queue<int> q;

	// Initialize visited table.
	for(int i = 0 ; i < vertices.size() ; i++)
		visited[i] = false;

	q.push(source);
	visited[source] = true;
	parentInfo[source] = -1;
	distanceInfo[source] = 0.0;

	while(!q.empty())
	{
		int currVertex = q.front();

		q.pop();

		for(int i = 0 ; i < vertices[currVertex]->outgoingEdges.size() ; i++)
		{
			int nextVertex = edges[vertices[currVertex]->outgoingEdges[i]]->v;

			if(!visited[nextVertex])
			{
				q.push(nextVertex);
				visited[nextVertex] = true;
				parentInfo[nextVertex] = currVertex;
				distanceInfo[nextVertex] = distanceInfo[currVertex] + edges[vertices[currVertex]->outgoingEdges[i]]->length;
			}
		}
	}

	for(int currVertex = dest ; currVertex != -1 ; currVertex = parentInfo[currVertex])
		path.push_back(currVertex);

	reverse(path.begin(), path.end());

	return path;
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
	Computes branch information vector, which holds segment connection information, abstracted away from 
	geometric details. 
*/
void Skeleton::computeBranchInfo(void)
{
	branchInfo.resize(vertices.size());

	for(int i = 0 ; i < vertices.size() ; i++)
	{
		if(vertices[i]->isRegular())
			continue;

		for(int j = 0 ; j < vertices[i]->outgoingEdges.size() ; j++)
		{
			int source = edges[vertices[i]->outgoingEdges[j]]->u;
			int dest = edges[vertices[i]->outgoingEdges[j]]->v;

			while(vertices[dest]->isRegular())
			{
				int tmp = dest;

				if(edges[vertices[dest]->outgoingEdges[0]]->v == source)
					dest = edges[vertices[dest]->outgoingEdges[1]]->v;
				else
					dest = edges[vertices[dest]->outgoingEdges[0]]->v;

				source = tmp;
			}

			branchInfo[i].push_back(dest);
		}
	}
}

/*
	Computes the segments of the skeleton in a list of list. Inner list contains vertex ids of vertexes in the piece in order. 
*/
void Skeleton::computeSegments(void)
{
	bool visitedEdges[edges.size()];

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
}

/*
	Computes distance between the skeleton and a point.
*/
double Skeleton::computeDistance2Skeleton(const Vector3f & point) const
{
	double d_min = FLT_MAX;

	for(int i = 0 ; i < edges.size() ; i += 2)
	{
		Vector3f u = vertices[edges[i]->u]->coords;
		Vector3f v = vertices[edges[i]->v]->coords;
		double d_curr = computeDistance2LineSeg(u, v, point);

		if(d_curr < d_min)
			d_min = d_curr;
	}

	return d_min;
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
	Computes average distance between this skeleton and skel. 
	TODO: How about squared distance? 
*/
double Skeleton::computeAvgDistance(double stepSize, const Skeleton *skel) const
{
	double avgDist = 0.0;
	unsigned long numSamples = 0;

	for(int i = 0 ; i < this->segments.size() ; i++)
	{
		tuple<double, int> currDist = this->computeDistInfoOfSegment(stepSize, this->segments[i], skel);

		avgDist += get<0>(currDist);
		numSamples += get<1>(currDist);
	}

	return avgDist / numSamples;
}

/*
	Samples the segment by stepSize. Evaluates distances of each sample to skel, sums them up and returns it as the first return value. 
	Second return value is the number of sample points. 
*/
tuple<double, int> Skeleton::computeDistInfoOfSegment(double stepSize, vector<int> segment, const Skeleton * & skel) const
{
	double totalDist = 0.0;
	Vector3f p_t;
	vector<double> distances;

	p_t = vertices[segment[0]]->coords;
	distances.push_back(skel->computeDistance2Skeleton(p_t));

	int i = 1;	// Index of the next vertex on the segment. 
	while(true)
	{
		// Move p_t on the arc as step size.
		double dist2Next = (vertices[segment[i]]->coords - p_t).norm();

		if(stepSize < dist2Next)
		{
			Vector3f dir = (vertices[segment[i]]->coords - p_t).normalized();

			p_t = p_t + dir * stepSize;
			distances.push_back(skel->computeDistance2Skeleton(p_t));
		}
		else if(stepSize == dist2Next)
		{
			p_t = vertices[segment[i]]->coords;
			distances.push_back(skel->computeDistance2Skeleton(p_t));
			i++;

			if(i == segment.size())
				break;
		}
		else	// stepSize >= dist2Next
		{
			p_t = vertices[segment[i]]->coords;
			i++;

			if(i == segment.size())
			{
				distances.push_back(skel->computeDistance2Skeleton(p_t));

				break;
			}

			double residual = stepSize - dist2Next;
			Vector3f dir = (vertices[segment[i]]->coords - p_t).normalized();

			p_t = p_t + dir * residual;
			distances.push_back(skel->computeDistance2Skeleton(p_t));
		}
	}

	// Integration over distances vector.
	for(i = 0 ; i < distances.size() ; i++)
		totalDist += distances[i] * distances[i];

	return make_tuple(totalDist, distances.size());
}

/*
	Computes Sup Int terms in Hausdorff distance with stepSize.
*/
double Skeleton::supInt(double stepSize, const Skeleton *skel) const
{
	double currMax;

	currMax = 0.0;
	for(int i = 0 ; i < this->segments.size() ; i++)
	{
		vector<Vector3f> samples = sampleSegment(stepSize, this->segments[i]);

		for(int j = 0 ; j < samples.size() ; j++)
		{
			double minDist = skel->computeDistance2Skeleton(samples[j]);

			if(minDist > currMax)
				currMax = minDist;
		}
	}

	return currMax;
}

/*
	Samples a skeleton segment (which is an arc) with a spacing of stepSize.
*/
vector<Vector3f> Skeleton::sampleSegment(double stepSize, vector<int> segment) const
{

	Vector3f p_t;
	vector<Vector3f> samples;

	p_t = vertices[segment[0]]->coords;
	samples.push_back(p_t);

	int i = 1;	// Index of the next vertex on the segment. 
	while(true)
	{
		// Move p_t on the arc as step size.
		double dist2Next = (vertices[segment[i]]->coords - p_t).norm();

		if(stepSize < dist2Next)
		{
			Vector3f dir = (vertices[segment[i]]->coords - p_t).normalized();

			p_t = p_t + dir * stepSize;
			samples.push_back(p_t);
		}
		else if(stepSize == dist2Next)
		{
			p_t = vertices[segment[i]]->coords;
			samples.push_back(p_t);
			i++;

			if(i == segment.size())
				break;
		}
		else	// stepSize >= dist2Next
		{
			p_t = vertices[segment[i]]->coords;
			i++;

			if(i == segment.size())
			{
				samples.push_back(p_t);

				break;
			}

			double residual = stepSize - dist2Next;
			Vector3f dir = (vertices[segment[i]]->coords - p_t).normalized();

			p_t = p_t + dir * residual;
			samples.push_back(p_t);
		}
	}

	return samples;
}
