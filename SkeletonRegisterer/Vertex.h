#ifndef _VERTEX_H_
#define _VERTEX_H_

#include <Eigen/Dense>
#include <vector>

using namespace Eigen;
using namespace std;	

class Vertex
{
public:
	int id;	// id of the vertex
	Vector3f coords;	// Coordinates of the vertex
	vector<int> outgoingEdges;	//List of outgoing edges

	Vertex(int id, Vector3f coords);	// Constructor
	int getDegree(void);	// Returns degree of the vertex
	bool isTerminal(void); // Checks if the vertex is a terminal vertex
	bool isRegular(void);	// Checks if the vertex is a regular vertex
	bool isJoint(void);	// Checks if the vertex is a joint vertex
};

#endif
