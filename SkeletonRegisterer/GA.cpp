#include "GA.h"

/*
	The constructor. 
	skel1: Curve skeleton (High-resolution)
	skel2: IK-skeleton (Low-resolution)
*/
GA::GA(const Skeleton * & skel1, const Skeleton * & skel2)
	: skel1(skel1), skel2(skel2)
{
	for(int i = 0 ; i < skel2->vertices.size() ; i++)
		if(skel2->vertices[i]->isTerminal())
			numTerminal++;
		else if(skel2->vertices[i]->isRegular())
			numRegular++;
		else
			numJoint++;

	for(int i = 0 ; i < skel1->vertices.size() ; i++)
		if(skel1->vertices[i]->isTerminal())
			terminalVertexPool.push_back(skel1->vertices[i]->id);
		else if(skel1->vertices[i]->isJoint())
			jointVertexPool.push_back(skel1->vertices[i]->id);

	genomeSize = skel2->vertices.size();

	population.resize(popSize);
	for(int i = 0 ; i  < popSize ; i++)
		population[i].genome.resize(genomeSize);

	totalFitness = 0.0;
}

/*
	Run the optimizer. Returns the optimum vertex matching between the IK-Skeleton and the curve-skeleton.
*/
Skeleton GA::run(void)
{
	printf("Initializing the population...\n");
	/*do
	{
		initializePopulation();	// Initialize the population.
	}
	while(totalFitness == 0.0);	// Never initialize a population in which fitness of all choromosomes are 0.0.
	printf("Population initialized.\n");*/

	/*fittestIdx = computeFittestIdx();
	fittestAllTime.push_back(population[fittestIdx]);*/
	fittestAllTime.push_back(population[0]);

	/*for(int i = 0 ; i < numIter ; i++)	// TODO: More advanced stop condition will be found later. 
	{
		// Print the current best fitness value reached. 
		printf("%d %lf %lf\n", i, fittestAllTime[0].fitness, computeAvgFitness());

		// First choose two parents to mate with fitness proportionate distribution.
		pair<int, int> parentIdxs = computeParents();
		int parentIdx1 = parentIdxs.first;
		int parentIdx2 = parentIdxs.second;

		// Do crossover operation.
		pair<Chromosome, Chromosome> offsprings = crossover(population[parentIdx1], population[parentIdx2]);
		Chromosome offspring1 = offsprings.first;
		Chromosome offspring2 = offsprings.second;

		// Do mutation operation.
		offspring1 = mutate(offspring1);
		offspring2 = mutate(offspring2);

		// Compute fitness values of new chromosomes. 
		offspring1.fitness = computeFitness(offspring1.genome);
		offspring2.fitness = computeFitness(offspring2.genome);

		// Choose two victims with inverse fitness proportionate distribution.
		pair<int, int> victimIdxs = computeVictims();
		int victimIdx1 = victimIdxs.first;
		int victimIdx2 = victimIdxs.second;

		// Update fittest chromosome of all time (if necessary)
		if(offspring1.fitness > fittestAllTime[0].fitness)
		{
			fittestAllTime.clear();
			fittestAllTime.push_back(offspring1);
		}
		else if(offspring1.fitness == fittestAllTime[0].fitness && isNewToFittestAllTime(offspring1.genome))
			fittestAllTime.push_back(offspring1);

		if(offspring2.fitness > fittestAllTime[0].fitness)
		{
			fittestAllTime.clear();
			fittestAllTime.push_back(offspring2);
		}
		else if(offspring2.fitness == fittestAllTime[0].fitness && isNewToFittestAllTime(offspring2.genome))
			fittestAllTime.push_back(offspring2);

		// Update total fitness of the population after removal of victim1.
		totalFitness -= population[victimIdx1].fitness;
		population[victimIdx1] = offspring1;
		totalFitness += population[victimIdx1].fitness;

		// Update index of the fittest chromosome after removal of victim1.
		if(victimIdx1 == fittestIdx)
			fittestIdx = computeFittestIdx();
		else if(offspring1.fitness > population[fittestIdx].fitness)
			fittestIdx = victimIdx1;

		// Update total fitness of the population after removal of victim2. 
		totalFitness -= population[victimIdx2].fitness;
		population[victimIdx2] = offspring2;
		totalFitness += population[victimIdx2].fitness;

		// Update intex of the fittest chromosome after removal of victim2.
		if(victimIdx2 == fittestIdx)
			fittestIdx = computeFittestIdx();
		else if(offspring2.fitness > population[fittestIdx].fitness)
			fittestIdx = victimIdx2;
	}

	// Choose the one with the smallest twist. 
	printf("Minimum twisted embedding is being computed...\n");
	int fittestWithMinTwist = computeMinimumTwistedSkeleton();
	printf("Minimum twisted embedding computed...\n");*/
	int fittestWithMinTwist = 0;


	/*printf("Number of optimum embeddings: %lu\n", fittestAllTime.size());

	for(int i = 0 ; i < fittestAllTime.size() ; i++)
	{
		fittestAllTime[i].fitness = computeFitness(fittestAllTime[i].genome);

		printf("%lf\n", fittestAllTime[i].fitness);
		for(int j = 0 ; j < genomeSize ; j++)
			printf("%d -> %d\n", j, fittestAllTime[i].genome[j]);

		printf("******\n");
	}*/

	return genomeToSkeleton(fittestAllTime[fittestWithMinTwist].genome);
}

/*
	Initializes the population. For each terminal and joint vertices in the IK skeleton, 
	it matches randomly terminal and joint vertices in the curve-skeleton. Terminal 
	vertices are only matched with terminal vertices and joint vertices are matched with 
	only joint vertices. 
*/
void GA::initializePopulation(void)
{
	mt19937 gen(rd());
	uniform_int_distribution<> disTerminal(0, terminalVertexPool.size() - 1);
	uniform_int_distribution<> disJoint(0, jointVertexPool.size() - 1);

	for(int i = 0 ; i < popSize ; i++)
	{
		bool usedTerminals[terminalVertexPool.size()];	// Do not match a terminal vertex in the curve-skeleton twice.
		bool usedJoints[jointVertexPool.size()];	// Do not match a joint vertex in the curve-skeleton twice.

		// Initialize the used vertex tables.
		for(int j = 0 ; j < terminalVertexPool.size() ; j++)
			usedTerminals[j] = false;
		for(int j = 0 ; j < jointVertexPool.size() ; j++)
			usedJoints[j] = false;

		// Fill the genome. 
		for(int j = 0 ; j < genomeSize ; j++)
		{
			if(skel2->vertices[j]->isTerminal())	// Match terminal vertex of IK-skeleton with a terminal vertex of curve-skeleton.
			{
				while(true)
				{
					int candidate = disTerminal(gen);

					if(!usedTerminals[candidate])	// Check if we used the candidate already. 
					{
						population[i].genome[j] = terminalVertexPool[candidate];
						usedTerminals[candidate] = true;

						break;
					}
				}
			}
			else if(skel2->vertices[j]->isRegular())	// Do not match regular vertices.
			{
				population[i].genome[j] = -1;
			}
			else	// Match joint vertex of IK-skeleton with a joint vertex of curve-skeleton.
			{
				while(true)
				{
					int candidate = disJoint(gen);

					if(!usedJoints[candidate])	// Check if we used the candidate already. 
					{
						population[i].genome[j] = jointVertexPool[candidate];
						usedJoints[candidate] = true;

						break;
					}
				}
			}
		}

		// Initialize the fitness of each chomosome too.
		population[i].fitness = computeFitness(population[i].genome);

		// Initialize the total fitness of the population. 
		totalFitness += population[i].fitness;
	}
}

/*
	Fitness proportionate selection.
*/
pair<int, int> GA::computeParents(void)
{
	int parent1 = 0;
	int parent2 = 0;
	double rand;
	mt19937 gen(rd());

	if(totalFitness == 0.0)
	{
		uniform_int_distribution<> disInt(0, popSize - 1); 

		do
		{
			parent1 = disInt(gen);
			parent2 = disInt(gen);
		}
		while(parent1 == parent2);

		return make_pair(parent1, parent2);
	}

	uniform_real_distribution<> dis(0.0, totalFitness);

	rand = dis(gen);

	for(double currSum = population[parent1].fitness ; currSum < rand ; currSum += population[parent1].fitness)
		parent1++;

	if(totalFitness - population[parent1].fitness == 0.0)
	{
		uniform_int_distribution<> disInt(0, popSize - 1);

		while(true)
		{
			parent2 = disInt(gen);

			if(parent1 != parent2)
				return make_pair(parent1, parent2);
		}
	}

	// May go to infinite loop if only one genome has a very large value than others. counter prevents that from happening. 
	int counter = 0;

	while(true)
	{
		parent2 = 0;

		rand = dis(gen);

		for(double currSum = population[parent2].fitness ; currSum < rand ; currSum += population[parent2].fitness)
			parent2++;

		if(parent1 != parent2)
			return make_pair(parent1, parent2);

		counter++;
		if(counter == 10)
		{
			uniform_int_distribution<> disInt(0, popSize - 1); 

			while(true)
			{
				parent2 = disInt(gen);

				if(parent1 != parent2)
					return make_pair(parent1, parent2);
			}
		}
	}
}

/*
	Inverse fitness proportionate selection.
*/
pair<int, int> GA::computeVictims(void)
{
	int victim1 = 0;
	int victim2 = 0;
	double invTotalFitness = 0.0;
	double invFitnessTable[popSize];
	double rand;
	vector<int> zeroFitnessChroms;	// Indices of chromosomes with zero fitness value. They should be victim directly. 
	mt19937 gen(rd());

	for(int i = 0 ; i < popSize ; i++)
	{
		if(population[i].fitness == 0.0)
		{
			zeroFitnessChroms.push_back(i);
			invFitnessTable[i] = INFINITY;
		}
		else
		{
			invFitnessTable[i] = 1.0 / population[i].fitness;
			invTotalFitness += invFitnessTable[i];
		}
	}

	if(zeroFitnessChroms.size() > 0)
	{
		uniform_int_distribution<> disInt(0, zeroFitnessChroms.size() - 1);

		victim1 = disInt(gen);

		if(zeroFitnessChroms.size() > 1)
		{
			while(true)
			{
				victim2 = disInt(gen);

				if(victim1 != victim2)
					return make_pair(victim1, victim2);
			}
		}
	}

	uniform_real_distribution<> dis(0.0, invTotalFitness);

	rand = dis(gen);

	if(zeroFitnessChroms.size() == 1)
	{
		double currSum;

		if(victim1 != 0)
		{
			victim2 = 0;
			currSum = invFitnessTable[0];
		}
		else
		{
			victim2 = 1;
			currSum = invFitnessTable[1];
		}

		while(currSum < rand)
		{
			victim2++;

			if(invFitnessTable[victim2] == INFINITY)
				continue;
			currSum += invFitnessTable[victim2];
		}

		return make_pair(victim1, victim2);
	}

	for(double currSum = invFitnessTable[victim1] ; currSum < rand ; currSum += invFitnessTable[victim1])
		victim1++;

	int counter = 0;

	while(true)
	{
		victim2 = 0;
		rand = dis(gen);

		for(double currSum = invFitnessTable[victim2] ; currSum < rand ; currSum += invFitnessTable[victim2])
			victim2++;

		if(victim1 != victim2)
			return make_pair(victim1, victim2);;

		// May go to infinite loop, counter prevents that from happening.
		counter++;
		if(counter == 10)
		{
			uniform_int_distribution<> disInt(0, popSize - 1); 

			while(true)
			{
				victim2 = disInt(gen);

				if(victim1 != victim2)
					return make_pair(victim1, victim2);
			}
		}
	}
}

/*
	Crossover operation will be done on the terminal vertices and joint vertices in an isolated way. 
	In other words, terminal vertices will only be matched with terminal vertices again, and joint 
	vertices will only be matched with joint vertices again, at the end of the operation. 
	Implements PMX (Partially-Mapped Crossover). 
*/
pair<GA::Chromosome, GA::Chromosome> GA::crossover(Chromosome chrom1, Chromosome chrom2)
{
	mt19937 gen(rd());
	uniform_int_distribution<> terminalCrossoverIndDist(1, numTerminal - 1);
	uniform_int_distribution<> jointCrossoverIndDist(1, numJoint - 1);
	uniform_real_distribution<> dis(0.0, 1.0);

	if(dis(gen) < crossoverProb)	// Crossover will be done with probability crossoverProb.
	{
		int numTerminalsSwapped = 0;	// Number of terminal genes swapped so far.  
		int numJointsSwapped = 0;	// Number of joint genes swapped so far.  

		int terminalXPoint = terminalCrossoverIndDist(gen);	// Choose the "cut point" for terminal genes randomly.
		int jointXPoint = jointCrossoverIndDist(gen);	// Choose the "cut point" for joint genes randomly. 

		Chromosome offspring1 = chrom1;
		Chromosome offspring2 = chrom2;

		// Not all indices will be related to a terminal or joint vertex. We should work only on those in an isolated manner. 
		for(int i = 0 ; i < genomeSize ; i++)
		{
			int duplicateInd;

			// Check if any interchange will occur.
			if((skel2->vertices[i]->isTerminal() && numTerminalsSwapped < terminalXPoint) || (skel2->vertices[i]->isJoint() && numJointsSwapped < jointXPoint))
			{
				// Write the gene from chrom2 to offspring1.
				duplicateInd = findGene(offspring1, chrom2.genome[i]);
				offspring1.genome[duplicateInd] = offspring1.genome[i];
				offspring1.genome[i] = chrom2.genome[i];

				// Write it from chrom1 to offspring2.
				duplicateInd = findGene(offspring2, chrom1.genome[i]);
				offspring2.genome[duplicateInd] = offspring2.genome[i];
				offspring2.genome[i] = chrom1.genome[i];

				if(skel2->vertices[i]->isTerminal())	// Check if a terminal gene has been swapped.
					numTerminalsSwapped++;
				else	// Joint gene has been swapped. 
					numJointsSwapped++;
			}
		}

		return make_pair(offspring1, offspring2);
	}

	return make_pair(chrom1, chrom2);
}

/*
	Mutation operation swaps randomly selected two genes. A terminal vertex is swapped only with a terminal vertex. 
	A joint vertex is swapped only with a joint vertex. Only one gene is swapped in a mutation operation. 
*/
GA::Chromosome GA::mutate(Chromosome chrom)
{
	mt19937 gen(rd());
	uniform_int_distribution<> mutationIndDist(0, genomeSize - 1);
	uniform_real_distribution<> dis(0.0, 1.0);

	if(dis(gen) < mutationProb)	// Mutation will be done with probability mutationProb. 
	{
		int mutInd1;	// Location of mutation will take place.
		int mutInd2;	// The other location of mutation.
		Chromosome mutant = chrom;	// The mutant.

		// Choose the locations to occur randomly.
		while(true)
		{
			mutInd1 = mutationIndDist(gen);
			mutInd2 = mutationIndDist(gen);

			if(mutInd1 == mutInd2)	// The mutation cannot occur at the same point.
				continue;

			if(skel2->vertices[mutInd1]->isTerminal() && skel2->vertices[mutInd2]->isTerminal())	// Terminal vertices can swap on their own.
				break;
			else if(skel2->vertices[mutInd1]->isJoint() && skel2->vertices[mutInd2]->isJoint())	// Joint vertices can swap on their own.
				break;
			else
				continue;
		}

		// Do mutation operation.
		mutant.genome[mutInd1] = chrom.genome[mutInd2];
		mutant.genome[mutInd2] = chrom.genome[mutInd1];

		return mutant;
	}

	return chrom;
}

/*
	Computes fitness value of a genome. 
*/
double GA::computeFitness(const vector<int> & genome) const
{
	double disSimTerm = 0.0;	// The lower the better. 
	double lengthDiffTerm = 0.0;	// The higher the better. 
	double changeAtJointAnglesTerm = 0.0;	// The lower the better. 
	double totalTwistTerm = 0.0;

	Skeleton tempSkel = genomeToSkeleton(genome);	// Embedded skeleton.

	// Hausdorff distance. Not used for now. 
	// fitness = 1.0 / skel1->computeHausdorffDistance(tempSkel);

	// Compute dissimilarity between the curve-skeleton and the embedded skeleton. 
	double disSim = skel1->computeDissimilarity(&tempSkel);

	// Compute length difference between the curve-skeleton and the embedded skeleton. 
	//double lengthDiff = skel1->totalLength - tempSkel.totalLength;

	// Compute change in pairwise angles at the joints considering the IK-skeleton and the embedded skeleton.
	//double changeAtJointAngles = computePairwiseAngleChange(&tempSkel);

	// Compute twist along skeleton segments with joints at both end considering the IK-skeleton and the embedded skeleton. 
	//double totalTwist = computeTotalTwist(&tempSkel);

	disSimTerm = 1.0 / disSim;
	//lengthDiffTerm = 1.0 * (lengthDiff <= 0.0 ? 0.0 : lengthDiff);
	//changeAtJointAnglesTerm = 1.0 / (changeAtJointAngles == 0.0 ? 0.1 : changeAtJointAngles);

	// printf("%lf %lf %lf\n", disSimTerm, lengthDiffTerm, changeAtJointAnglesTerm);

	return disSimTerm + lengthDiffTerm + changeAtJointAnglesTerm + totalTwistTerm;
}

/*
	Returns index of the fittest chromosome in the population. 
*/
int GA::computeFittestIdx(void) const
{
	double currMax = -1.0;
	int currFittest = -1;

	for(int i = 0 ; i < popSize ; i++)
		if(population[i].fitness > currMax)
		{
			currMax = population[i].fitness;
			currFittest = i;
		}

	return currFittest;
}

/*
	Computes average fitness of the population. 
*/
double GA::computeAvgFitness(void) const
{
	return totalFitness / popSize;
}

/*
	Returns position of gene inside chrom.genome vector.
*/
int GA::findGene(Chromosome & chrom, int gene)
{
	for(int i = 0 ; i < genomeSize ; i++)
		if(gene == chrom.genome[i])
			return i;

	return -1;
}

/*
	Checks whether the genome exists in fittestAllTime or not. 
*/
bool GA::isNewToFittestAllTime(const vector<int> & genome) const
{
	for(int i = 0 ; i < fittestAllTime.size() ; i++)
		if(genome == fittestAllTime[i].genome)
			return false;

		return true;	
}

/*
	Takes a genome as an argument, which supplies terminal and joint vertex matching 
	between the IK skeleton and curve-skeleton. Returns the skeleton it represents 
	by computing regular vertex positions. 
*/
Skeleton GA::genomeToSkeleton(const vector<int> & genome) const
{
	Skeleton representedSkel;

	// Branch info is always same with the IK-skeleton.
	representedSkel.branchInfo = skel2->branchInfo;

	// Segments are always same with the IK-skeleton.
	representedSkel.segments = skel2->segments;

	// Compute edge information of the skeleton. Edge lengths will be computed after optimum positions of vertices will be computed. 
	for(int i = 0 ; i < skel2->edges.size() ; i++)
		representedSkel.edges.push_back(new Edge(skel2->edges[i]->id, skel2->edges[i]->u, skel2->edges[i]->v, skel2->edges[i]->opp, -1.0));	

	for(int i = 0 ; i < genome.size() ; i++)
	{
		if(genome[i] == -1)	// Placing regular vertex. 
		{
			tuple<int, double, int, double> segmentInfo = computeSegmentEnds(i, skel2);

			int end1 = genome[get<0>(segmentInfo)];	// We are working with skel1 here. 
			int end2 = genome[get<2>(segmentInfo)];

			double geo1 = get<1>(segmentInfo);	// The geodesic distance between i and end1. 
			double geo2 = get<3>(segmentInfo);	// The geodesic distance between i and end2.
			double ratio = geo1 / (geo1 + geo2);	// The ratio should be kept. 

			Vector3f optPosition = placeRegularVertex(end1, end2, ratio);

			representedSkel.vertices.push_back(new Vertex(i, optPosition));
		}
		else	// Placing non-regular vertex. 
			representedSkel.vertices.push_back(new Vertex(i, skel1->vertices[genome[i]]->coords));
	}

	int shortestEdgeIdx = -1;
	double minEdgeLength = FLT_MAX;
	for(int i = 0 ; i < representedSkel.edges.size() ; i += 2)
	{
		int v1 = representedSkel.edges[i]->u;
		int v2 = representedSkel.edges[i]->v;
		double edgeLength = (representedSkel.vertices[v1]->coords - representedSkel.vertices[v2]->coords).norm();

		representedSkel.totalLength += edgeLength;
		representedSkel.edges[i]->length = edgeLength;
		representedSkel.edges[i + 1]->length = edgeLength;

		if(edgeLength < minEdgeLength)
		{
			minEdgeLength = edgeLength;
			shortestEdgeIdx = i;
		}
	}
	representedSkel.shortestEdgeIdx = shortestEdgeIdx;

	// Copy the outgoing edge information as is from the IK-skeleton.
	for(int i = 0 ; i < representedSkel.vertices.size() ; i++)
		representedSkel.vertices[i]->outgoingEdges = skel2->vertices[i]->outgoingEdges;

	return representedSkel;
}

/*
	Takes vertexId as argument. Computes ids of end points of vertices at the two ends of the segment 
	on which the vertex with vertexId lies. 
	Returns, id of end vertex and its geodesic distance to vertexId, id for the other end vertex and 
	its geodesic distance to vertexId. 
*/
tuple<int, double, int, double> GA::computeSegmentEnds(int vertexId, const Skeleton * skel) const
{	
	int prevEndVertex1 = vertexId;
	int endVertex1 = skel->edges[skel->vertices[prevEndVertex1]->outgoingEdges[0]]->v;
	double geo1 = skel->edges[skel->vertices[prevEndVertex1]->outgoingEdges[0]]->length;

	while(skel->vertices[endVertex1]->isRegular())
	{
		int tmp = endVertex1;

		if(skel->edges[skel->vertices[endVertex1]->outgoingEdges[0]]->v == prevEndVertex1)
		{
			geo1 += skel->edges[skel->vertices[endVertex1]->outgoingEdges[1]]->length;
			endVertex1 = skel->edges[skel->vertices[endVertex1]->outgoingEdges[1]]->v;
		}
		else
		{
			geo1 += skel->edges[skel->vertices[endVertex1]->outgoingEdges[0]]->length;
			endVertex1 = skel->edges[skel->vertices[endVertex1]->outgoingEdges[0]]->v;
		}

		prevEndVertex1 = tmp;
	}

	int prevEndVertex2 = vertexId;
	int endVertex2 = skel->edges[skel->vertices[prevEndVertex2]->outgoingEdges[1]]->v;
	double geo2 = skel->edges[skel->vertices[prevEndVertex2]->outgoingEdges[1]]->length;
	while(skel->vertices[endVertex2]->isRegular())

	{
		int tmp = endVertex2;

		if(skel->edges[skel->vertices[endVertex2]->outgoingEdges[0]]->v == prevEndVertex2)
		{
			geo2 += skel->edges[skel->vertices[endVertex2]->outgoingEdges[1]]->length;
			endVertex2 = skel->edges[skel->vertices[endVertex2]->outgoingEdges[1]]->v;
		}
		else
		{
			geo2 += skel->edges[skel->vertices[endVertex2]->outgoingEdges[0]]->length;
			endVertex2 = skel->edges[skel->vertices[endVertex2]->outgoingEdges[0]]->v;
		}

		prevEndVertex2 = tmp;
	}

	return make_tuple(endVertex1, geo1, endVertex2, geo2);
}

/*
	Computes position of the regular vertex on the skeleton segment starting with terminal vertex end1 and 
	ending with terminal vertex end2. end1 corresponds to 0.0 and end2 corresponds to 1.0 in parametric segment 
	representation. Returns the position when we move from end1 to end2 on the skeleton segment with t=ratio. 
*/
Vector3f GA::placeRegularVertex(int end1, int end2, double ratio) const
{
	int v1;
	int v2;
	double distanceInfo[skel1->vertices.size()];
	Vector3f optPos;

	vector<int> path = skel1->computeShortestPath_BFS(end1, end2, distanceInfo);	// Compute shortest path between end1 and end2.

	double totalDistance = distanceInfo[end2];	// Total distance is the distance between end1 and end2.

	for(int i = 0 ; i < path.size() - 1 ; i++)
	{
		v1 = path[i];
		v2 = path[i + 1];

		double t = distanceInfo[v2] / totalDistance;	// t is the normalized distance of v2.
		if(t >= ratio)
			break;
	}

	if(v1 >= skel1->vertices.size())
		printf("Segment!!!\n");
	double t1 = distanceInfo[v1] / totalDistance;
	double t2 = distanceInfo[v2] / totalDistance;
	double l1 = distanceInfo[v1];
	double l2 = distanceInfo[v2];
	double residual = (ratio - t1) / (t2 - t1);

	optPos = skel1->vertices[v1]->coords + (skel1->vertices[v2]->coords - skel1->vertices[v1]->coords) * residual * (l2 - l1);

	return optPos;
}

/*
	Computes change in pairwise angles at the joints considering the IK-skeleton and skel.
*/
double GA::computePairwiseAngleChange(const Skeleton * skel) const
{
	double totalChange = 0.0;

	for(int i = 0 ; i < skel->vertices.size() ; i++)
		if(skel->vertices[i]->isJoint())
			totalChange += pairwiseAngleChangePerJoint(skel, i);

	return totalChange;
}

/*
	Computes change in pairwise angles in the joint vertex idx considering the IK-skeleton and skel (The embedded skeleton). 
*/
double GA::pairwiseAngleChangePerJoint(const Skeleton * skel, int idx) const
{
	double changeAtJoint = 0.0;
	Vector3f p_jnt_IK = skel2->vertices[idx]->coords;	// Joint vertex position in the IK-skeleton.
	Vector3f p_jnt_emb = skel->vertices[idx]->coords;	// Joint vertex position in the embedded skeleton.

	for(int i = 0 ; i < skel->vertices[idx]->outgoingEdges.size() - 1 ; i++)
	{
		// Position of the vertex at the ith IK-skeleton segment.
		Vector3f p_seg_IK_1 = skel2->vertices[skel2->edges[skel2->vertices[idx]->outgoingEdges[i]]->v]->coords;
		// Position of the vertex at the ith embedded segment.
		Vector3f p_seg_emb_1 = skel->vertices[skel->edges[skel->vertices[idx]->outgoingEdges[i]]->v]->coords;

		// Normalized vector associated with the ith outgoing edge in the IK-skeleton. 
		Vector3f e_IK_1 = (p_seg_IK_1 - p_jnt_IK).normalized();
		// Normalized vector associated with the ith outgoing edge in the embedded skeleton.
		Vector3f e_emb_1 = (p_seg_emb_1 - p_jnt_emb).normalized();

		for(int j = i + 1 ; j < skel->vertices[idx]->outgoingEdges.size() ; j++)
		{
			// Position of the vertex at the jth IK-skeleton segment.
			Vector3f p_seg_IK_2 = skel2->vertices[skel2->edges[skel2->vertices[idx]->outgoingEdges[j]]->v]->coords;
			// Position of the vertex at the jth embedded skeleton segment.
			Vector3f p_seg_emb_2 = skel->vertices[skel->edges[skel->vertices[idx]->outgoingEdges[j]]->v]->coords;

			// Normalized vector associated with the jth outgoing edge in the IK-skeleton. 
			Vector3f e_IK_2 = (p_seg_IK_2 - p_jnt_IK).normalized();
			// Normalized vector associated with the jth outgoing edge in the embedded skeleton. 
			Vector3f e_emb_2 = (p_seg_emb_2 - p_jnt_emb).normalized();

			double cosine_IK = e_IK_1.dot(e_IK_2);
			double cosine_em = e_emb_1.dot(e_emb_2);

			changeAtJoint += fabs(cosine_IK - cosine_em);
		}
	}

	return changeAtJoint;
}

/*
	Returns the index (in fittestAllTime) of the minimum twisted skeleton in the fittest genes. 
	TODO: Implement this one. 
*/
int GA::computeMinimumTwistedSkeleton(void)
{
	int minTwistedIdx = 0;

	return minTwistedIdx;
}

/*
	Compute twist along skeleton segments with joints at both end considering the IK-skeleton and skel.
	TODO: The segment having joint vertices at both ends may contain regular vertices. The method should also 
	handle this case.  
*/
double GA::computeTotalTwist(const Skeleton * skel) const 
{
	double totalTwist = 0.0;

	// Find segments having joint vertices at both ends. 
	for(int i = 0 ; i < skel2->segments.size() ; i++)
	{
		int vert1 = skel2->segments[i][0];
		int vert2 = skel2->segments[i][skel2->segments[i].size() - 1];

		if(skel2->vertices[vert1]->isJoint() && skel2->vertices[vert2]->isJoint())
			totalTwist += twistPerSegment(skel, vert1, vert2);
	}

	return totalTwist;
}

double GA::twistPerSegment(const Skeleton * skel, int vert1, int vert2) const
{
	int numEdgePairs = 0;
	double avgTwistInSegment = 0.0;
	Vector3f normal_IK = (skel2->vertices[vert1]->coords - skel2->vertices[vert2]->coords).normalized();
	Vector3f normal_emb = (skel->vertices[vert1]->coords - skel->vertices[vert2]->coords).normalized();

	for(int i = 0 ; i < skel2->vertices[vert1]->outgoingEdges.size() ; i++)
	{
		Vector3f s_i_IK = skel2->vertices[skel2->edges[skel2->vertices[vert1]->outgoingEdges[i]]->v]->coords - skel2->vertices[vert1]->coords;

		if(fabs(normal_IK.dot(s_i_IK.normalized())) == 1.0)
			continue;

		Vector3f s_i_emb = skel->vertices[skel->edges[skel->vertices[vert1]->outgoingEdges[i]]->v]->coords - skel->vertices[vert1]->coords;

		if(fabs(normal_emb.dot(s_i_emb.normalized())) == 1.0)
			continue;

		Vector3f s_i_IK_ = normal_IK.dot(s_i_IK) * normal_IK;
		Vector3f s_i_IK__ = (s_i_IK - s_i_IK_).normalized();

		Vector3f s_i_emb_ = normal_emb.dot(s_i_emb) * normal_emb;
		Vector3f s_i_emb__ = (s_i_emb - s_i_emb_).normalized();

		for(int j = 0 ; j < skel2->vertices[vert2]->outgoingEdges.size() ; j++)
		{
			Vector3f s_j_IK = skel2->vertices[skel2->edges[skel2->vertices[vert2]->outgoingEdges[j]]->v]->coords - skel2->vertices[vert2]->coords;

			if(fabs(normal_IK.dot(s_j_IK.normalized())) == 1.0)
				continue;

			Vector3f s_j_emb = skel->vertices[skel->edges[skel->vertices[vert2]->outgoingEdges[j]]->v]->coords - skel->vertices[vert2]->coords;

			if(fabs(normal_emb.dot(s_j_emb.normalized())) == 1.0)
				continue;

			Vector3f s_j_IK_ = normal_IK.dot(s_j_IK) * normal_IK;
			Vector3f s_j_IK__ = (s_j_IK - s_j_IK_).normalized();

			Vector3f s_j_emb_ = normal_emb.dot(s_j_emb) * normal_emb;
			Vector3f s_j_emb__ = (s_j_emb - s_j_emb_).normalized();

			avgTwistInSegment += fabs(s_i_IK__.dot(s_j_IK__) - s_i_emb__.dot(s_j_emb__));

			numEdgePairs++;
		}
	}

	return avgTwistInSegment / numEdgePairs;
}
