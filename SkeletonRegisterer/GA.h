#ifndef _GA_H_
#define _GA_H_

#include <cfloat>
#include <Eigen/Dense>
#include <random>
#include <vector>
#include <tuple>
#include <utility>

#include "Skeleton.h"
#include "Vertex.h"

using namespace Eigen;
using namespace std;

class GA
{
public: 
	GA(const Skeleton * & skel1, const Skeleton * & skel2);
	Skeleton run(void);

private:
	struct Chromosome
	{
		vector<int> genome;
		double fitness;
	};

	int popSize = 5;			// Population size
	int numTerminal = 0;		// Number of terminal vertices in skel2 (IK Skeleton)
	int numJoint = 0;			// Number of joint vertices in skel2 (IK Skeleton)
	int numRegular = 0;			// Number of regular vertice in skel2 (IK Skeleton)
	int genomeSize;				// Size of a genome
	int fittestIdx;				// Index of the fittest chromosome
	int numIter = 1;			// Number of iterations to stop
	float crossoverProb = 0.8;	// Probability of crossover occurrence
	float mutationProb = 0.4;	// Probability of mutation occurrence
	double totalFitness;		// Total fitness of the population
	vector<int> terminalVertexPool;	// Set of terminal vertices in skel1 (curve-skeleton)
	vector<int> jointVertexPool;	// Set of joint vertices in skel1 (curve-skeleton)
	vector<Chromosome> population;	// Main population
	vector<Chromosome> fittestAllTime;	// Fittest chromosomes encountered so far
	random_device rd;
	const Skeleton *skel1;	// High resolution
	const Skeleton *skel2;	// Low resolution

	void initializePopulation(void);
	pair<int, int> computeParents(void);
	pair<int, int> computeVictims(void);
	pair<Chromosome, Chromosome> crossover(Chromosome chrom1, Chromosome chrom2);
	Chromosome mutate(Chromosome chrom);
	double computeFitness(const vector<int> & genome) const;	// Computes fitness of a genome
	int computeFittestIdx(void) const;	// Computes the fittest index of the fittest chromosome in the population
	double computeAvgFitness(void) const;	// Computes average fitness of the population. 
	int findGene(Chromosome & chrom, int gene);	// Returns index of gene in chrom 
	bool isNewToFittestAllTime(const vector<int> & genome) const; // Checks whether the genome exists in fittestAllTime or not. 
	Skeleton genomeToSkeleton(const vector<int> & genome) const; // Constructs the skeleton that a genome represents.
	tuple<int, double, int, double> computeSegmentEnds(int vertexId, const Skeleton * skel) const;	// Computes ids of end vertices of the segment vertexId lies. 
	Vector3f placeRegularVertex(int end1, int end2, double ratio) const;	// Places regular vertex into segment according to ratio.
	double computePairwiseAngleChange(const Skeleton * skel) const; // Computes change in pairwise angles at the joints considering the IK-skeleton and skel.
	double pairwiseAngleChangePerJoint(const Skeleton * skel, int idx) const; // Computes change in pairwise angles in the joint vertex idx. 
	int computeMinimumTwistedSkeleton(void); // Returns the index (in fittestAllTime) of the minimum twisted skeleton in the fittest genes. 
	double computeTotalTwist(const Skeleton * skel) const; // Compute twist along skeleton segments with joints at both end considering the IK-skeleton and skel. 
	double twistPerSegment(const Skeleton * skel, int vert1, int vert2) const;
};

#endif
