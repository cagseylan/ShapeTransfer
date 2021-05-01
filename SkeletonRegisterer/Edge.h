#ifndef _EDGE_H_
#define _EDGE_H_

class Edge
{
public: 
	int id;	// id of the edge.
	int u;	// id of the source vertex.
	int v;	// id of the dest. vertex.
	int opp = -1;	// id of the opposite half-edge.
	double length;	// length of the edge. 

	Edge(int id, int u, int v, int opp, double length);
};

#endif
