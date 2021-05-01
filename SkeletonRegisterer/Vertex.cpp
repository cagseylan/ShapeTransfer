#include "Vertex.h"

Vertex::Vertex(int id, Vector3f coords) 
	: id(id), coords(coords)
{}

int Vertex::getDegree(void)
{
	return outgoingEdges.size();
}

bool Vertex::isTerminal(void)
{
	return outgoingEdges.size() == 1;
}

bool Vertex::isRegular(void)
{
	return outgoingEdges.size() == 2;
}

bool Vertex::isJoint(void)
{
	return outgoingEdges.size() > 2;
}
