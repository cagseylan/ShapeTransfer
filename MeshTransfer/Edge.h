#ifndef _EDGE_H_
#define _EDGE_H_

class Edge
{
public: 
	int id;	// id of the edge. 
	int u;	// id of the source vertex. 
	int v;	// id of the dest. vertex. 
	int opp = -1;	// id of the opposite half-edge. 
	int faceId; // Id of the triangle to which the edge corresponds. 
	float length; // Length of the edge. 
	float w = 0.0; // Weight to be used in W and L. 

	Edge(int id, int u, int v, int opp, int faceId);
};

#endif
