//	Author: Christopher Lee Jackson & Jason Jones
//	Course: CMPSC463
//  Problem: 3-1
//  Description:

#ifndef __ABDBMST_H_INCLUDED__  
#define __ABDBMST_H_INCLUDED__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <vector>
#include "Graph.h"
#include "BinaryHeap.h"
#include <cmath>
#include <cstring>

using namespace std;

bool asc_cmp_plevel(Edge *a, Edge *b);
bool des_cmp_cost(Edge *a, Edge *b);
bool asc_src(Edge *a, Edge *b);
void printEdge(Edge* e);

class Ant;
class Range;

class Ant {
public:
	int data; // initial vertex ant started on
	Vertex *location;
	vector<int> *visited;
};

class Range {
public:
	double low;
	double high;
	Edge *assocEdge;
};

class abdbmst {
private:
	//	
	static const double P_UPDATE_EVAP;
	static const double P_UPDATE_ENHA;
	static const int TABU_MODIFIER;
	static const int MAX_CYCLES;
	//	Class Variables
	double evap_factor, enha_factor, loopCount, maxCost, minCost, bestCost;
	int cycles, totalCycles, d;
	unsigned int MAX_TREE_SIZE;
	Graph* g;
	vector<Edge*> best;
	//	Prototypes
	void processEFile(char* file, int d);
	void processFileOld(char* file, int d);
	void processRFile(char* file, int d);
	vector<Edge*> AB_DBMST();
	vector<Edge*> treeConstruct();
	void move(Ant *a);
	void updatePheromonesPerEdge();
	void updatePheromonesGlobal(bool improved);

	void compute();	
public:
	//	Constructors
	abdbmst(Graph *gin, int din, double minin, double maxin);
	//	Prototypes

};

const double abdbmst::P_UPDATE_EVAP = 0.95;
const double abdbmst::P_UPDATE_ENHA = 1.05;
const int abdbmst::TABU_MODIFIER = 5;
const int abdbmst::MAX_CYCLES = 2500;

abdbmst::abdbmst(Graph *gin, int din, double minin, double maxin ) {
	cycles = 1;
	totalCycles = 1;
	minCost = minin;
	maxCost = maxin;
	evap_factor = 0.5;
	enha_factor = 1.5;
	loopCount = 0;
	g = gin;
	d = din;
	compute();
	MAX_TREE_SIZE = g->getCount() - 1;
	bestCost = std::numeric_limits<double>::infinity();
}

void abdbmst::compute() {
	cout << "Tree Size: " << g->getCount() << endl;
	cout << "Diameter Bound: " << d << endl;
	cout << "Best Size: " << best.size() << endl;
	cout << "Best Cost: " << bestCost << endl;
	cout << "running.....";
	AB_DBMST();
	cout << "done\n";
	cout << "Tree Size: " << g->getCount() << endl;
	cout << "Diameter Bound: " << d << endl;
	cout << "Best Size: " << best.size() << endl;
	cout << "Best Cost: " << bestCost << endl;
	//sort(best.begin(), best.end(), asc_src);
	cout << "Best Tree num edges: " << best.size() << endl;
	for_each(best.begin(), best.end(), printEdge);
	//best.clear();
}

/*
*
*
*	For Each, Sort - Helper Functions
*
*/
bool asc_cmp_plevel(Edge *a, Edge *b) {
	return (a->pLevel < b->pLevel);
}

bool des_cmp_cost(Edge *a, Edge *b) {
	return (a->weight > b->weight);
}

bool asc_src(Edge* a, Edge* b) {
	return (a->getSource(NULL)->data < b->getSource(NULL)->data);
}

void printEdge(Edge* e) {
	cout << e->getSource(NULL)->data << " " << e->getDestination(NULL)->data << " " << e->weight << " " << e->pLevel << endl;
}

vector<Edge*> abdbmst::AB_DBMST() {
	//	Local Declerations
	double treeCost = 0;
	const int s = 75;
	vector<Edge*> current;
	Vertex *vertWalkPtr;
	Edge *edgeWalkPtr;
	vector<Ant*> ants;
	vector<Edge*>::iterator e, ed;
	Ant *a;
	//	Assign one ant to each vertex
	vertWalkPtr = g->getFirst();
	for (unsigned int i = 0; i < g->getCount(); i++) {
		Ant *a = new Ant;
		a->data = i +1;
		a->location = vertWalkPtr;
		a->visited = new vector<int>(g->getCount(), 0);
		ants.push_back(a);
		//	Initialize pheremone level of each edge, and set pUdatesNeeded to zero
		for ( e = vertWalkPtr->edges.begin() ; e < vertWalkPtr->edges.end(); e++ ) {
			edgeWalkPtr = *e;
			if (edgeWalkPtr->getSource(NULL) == vertWalkPtr) {
				edgeWalkPtr->pUpdatesNeeded = 0;
				edgeWalkPtr->pLevel = (maxCost - edgeWalkPtr->weight) + ((maxCost - minCost) / 3);
			}
		}
		//	Done with this vertex's edges; move on to next vertex
		vertWalkPtr = vertWalkPtr->pNextVert;
	}
	while (totalCycles <= 10000 && cycles <= MAX_CYCLES) { 
		//cerr << "Cycle" << totalCycles << endl;
		if(totalCycles % 100 == 0) 
			cerr << "CYCLE " << totalCycles << endl;
		//	Exploration Stage
		for (int step = 1; step <= s; step++) {
			if (step == s/3 || step == (2*s)/3) {
				updatePheromonesPerEdge();
			}
			for (unsigned int j = 0; j < g->getCount(); j++) {
				a = ants[j];
				move(a);
			}
			if ( step % TABU_MODIFIER == 0 ) {
				for(unsigned int w = 0; w < g->getCount(); w++) {
					ants[w]->visited->assign(g->getCount(), 0); //  RESET VISITED FOR EACH ANT (TABU)
				}
			}
		}
		for(unsigned int w = 0; w < g->getCount(); w++) {
			ants[w]->visited->assign(g->getCount(), 0); //  RESET VISITED FOR EACH ANT
		}
		updatePheromonesPerEdge();
		//	Tree Construction Stage
		current = treeConstruct();
		//	Get new tree cost
		for ( ed = current.begin() ; ed < current.end(); ed++ ) {
			edgeWalkPtr = *ed;
			treeCost+=edgeWalkPtr->weight;
		}
		if (treeCost < bestCost) {
			cout << "FOUND NEW BEST" << " at cycle: " << totalCycles << endl;
			best = current;
			bestCost = treeCost;
			if (totalCycles != 1)
				cycles = 0;
		} 
		if (cycles % 100 == 0) {
			updatePheromonesGlobal(false);
		} else {
			updatePheromonesGlobal(true);
		}
		if (totalCycles % 500 == 0) {
			evap_factor *= P_UPDATE_EVAP; 
			enha_factor *= P_UPDATE_ENHA; 
		}
		totalCycles++;
		cycles++;
		treeCost = 0;
	}
	//cout << "Cycles: " << totalCycles << endl;
	cout << bestCost << endl;
	cycles =0;
	totalCycles = 0;
	return best;
}

void abdbmst::updatePheromonesGlobal(bool improved) {
	//	Local Variables
	srand((unsigned) time(NULL));
	double pMax = 1000*((maxCost - minCost) + (maxCost - minCost) / 3);
	double pMin = (maxCost - minCost)/3;
	Edge *e;
	double XMax = 0.3;
	double XMin = 0.1;
	double rand_evap_factor;
	double IP;
	vector<Edge*>::iterator ex;
	//	For each edge in the best tree update pheromone levels
	for ( ex = best.begin() ; ex < best.end(); ex++ ) {
		e = *ex;
		IP = (maxCost - e->weight) + ((maxCost - minCost) / 3);
		if (improved) {
			//	IMPROVEMENT so Apply Enhancement
			e->pLevel = enha_factor*e->pLevel;
		} else {
			//	NO IMPROVEMENTS so Apply Evaporation
			rand_evap_factor = XMin + rand() * (XMax - XMin) / RAND_MAX;
			e->pLevel = rand_evap_factor*e->pLevel;
		}
		//	Check if fell below minCost or went above maxCost
		if (e->pLevel > pMax) {
			e->pLevel = pMax - IP;
		} else if (e->pLevel < pMin) {
			e->pLevel = pMin + IP;
		}
	}
}

void abdbmst::updatePheromonesPerEdge() {
	//	Local Variables
	Vertex *vertWalkPtr = g->getFirst();
	double pMax = 1000*((maxCost - minCost) + (maxCost - minCost) / 3);
	double pMin = (maxCost - minCost)/3;
	double IP;
	vector<Edge*>::iterator ex;
	Edge *edgeWalkPtr;

	while (vertWalkPtr) {
		for ( ex = vertWalkPtr->edges.begin() ; ex < vertWalkPtr->edges.end(); ex++ ) {
			edgeWalkPtr = *ex;
			if (edgeWalkPtr->getSource(NULL) == vertWalkPtr) {
				IP = (maxCost - edgeWalkPtr->weight) + ((maxCost - minCost) / 3);
				edgeWalkPtr->pLevel = (1 - evap_factor)*(edgeWalkPtr->pLevel)+(edgeWalkPtr->pUpdatesNeeded * IP);
				if (edgeWalkPtr->pLevel > pMax) {
					edgeWalkPtr->pLevel = pMax - IP;
				} else if (edgeWalkPtr->pLevel < pMin) {
					edgeWalkPtr->pLevel = pMin + IP;
				}
				//	Done updating this edge reset multiplier
				edgeWalkPtr->pUpdatesNeeded = 0;
			}
		}
		vertWalkPtr = vertWalkPtr->pNextVert;
	}
}

vector<Edge*> abdbmst::treeConstruct() {
	//	Local Variables
	vector<Edge*> v, c, tree, possConn;
	const int HUBS_NEEDED = d - 1;
	vector<Hub*> hubs, treeHubs;
	Vertex *vertWalkPtr, *vert, *v1, *v2;
	vector<Hub*> possVerts;
	Hub *pHub, *h;
	Edge *pE, *pEdge, *edgeWalkPtr;
	int vertIndex;
	int numHubs = 0;
	unsigned int treeCount = 0;
	vector<Edge*>::iterator iedge1, iedge2, iedge3, ie;
	vector<Hub*>::iterator ihubs1, ihubs2;
	Hub *highHub = NULL;
	vector<int> uf( g->getCount()+1 , 0 );
	BinaryHeap* heap;
	//	Put all edges into a vector
	vertWalkPtr = g->getFirst();
	while (vertWalkPtr) {
		vertWalkPtr->treeDegree = 0;
		vertWalkPtr->inTree = false;
		vertWalkPtr->isConn = false;
		for ( ie = vertWalkPtr->edges.begin() ; ie < vertWalkPtr->edges.end(); ie++ ) {
			edgeWalkPtr = *ie;
			//	Dont want duplicate edges in listing
			if (edgeWalkPtr->getSource(NULL) == vertWalkPtr) {
				edgeWalkPtr->inTree = false;
				v.push_back(edgeWalkPtr);
			}
		}
		vertWalkPtr = vertWalkPtr->pNextVert;
	}
	//	Sort edges in ascending order based upon pheromone level
	sort(v.begin(), v.end(), asc_cmp_plevel);
	//	Select 5n edges from the end of v( the highest pheromones edges) and put them into c.
	for (unsigned int i = 0; i < 5*g->getCount(); i++) {
		if (v.empty()) {
			break;
		}
		c.push_back(v.back());
		v.pop_back();
	}
	//	Sort edges in descending order based upon cost
	sort(c.begin(), c.end(), des_cmp_cost);
	//  Fill vector of Hubs
	vert = g->getFirst();
	for(unsigned int index = 0; index < g->getCount(); index++) {
		hubs.push_back(new Hub());
		hubs[index]->vertId = index + 1;
		hubs[index]->vert = vert;
		vert = vert->pNextVert;
	}
	//  Now get d - 1 hubs
	while(numHubs < HUBS_NEEDED && treeCount != MAX_TREE_SIZE) {
		if(!c.empty()){
			//  Get Degree of each vertice in candidate set
			for(iedge1 = c.begin(); iedge1 < c.end(); iedge1++) {
				pEdge = *iedge1;
				//  Handle Source
				vertIndex = pEdge->getSource(NULL)->data; // the vertice number uniquely identifies each vertice
				hubs[vertIndex - 1]->edges.push_back(pEdge);
				//  Handle Destination
				vertIndex = pEdge->getDestination(NULL)->data; // the vertice number uniquely identifies each vertice
				hubs[vertIndex - 1]->edges.push_back(pEdge);
			}
			//  Put Potential hubs in to heap
			//  First get rid of vertices with zero edges from candidate set
			for(ihubs2 = hubs.begin(); ihubs2 < hubs.end(); ihubs2++) {
				pHub = *ihubs2;
				if(pHub->edges.size() != 0) {
					possVerts.push_back(pHub);
				}
			}
			heap = new BinaryHeap( possVerts );
			//  Get highest degree v to make our initial hub (should be top of heap)
			highHub = heap->deleteMax();
			numHubs++;
			treeHubs.push_back(highHub);
			//  Add all edges in highHub to tree
			for(iedge1 = highHub->edges.begin(); iedge1 < highHub->edges.end(); iedge1++) {
				pEdge = *iedge1;
				if(!pEdge->inTree) {
					pEdge->getDestination(NULL)->inTree = true;
					pEdge->getSource(NULL)->inTree = true;
					pEdge->inTree = true;
					tree.push_back(pEdge);
					treeCount++;
				}
			}
			//  Update potential connector edges
			if ( numHubs > 1 ) {
				for(int i = numHubs - 1; i >= 1; i--) {
					v1 = highHub->vert; 
					v2 = treeHubs[i]->vert;
					for(iedge3 = v2->edges.begin(); iedge3 < v2->edges.end(); iedge3++) {
						pEdge = *iedge3;
						if(pEdge->getDestination(NULL)->data ==  v1->data) {
							possConn.push_back(pEdge);
						}
					}
				}
			}
				//  Get rid of edges that are already in tree or that would cause a loop
			for(iedge1 = highHub->edges.begin(); iedge1 < highHub->edges.end(); iedge1++) {
				pEdge = *iedge1;
				//  Update Source Vertex
				h = hubs[pEdge->getSource(NULL)->data - 1];
				for(iedge2 = h->edges.begin() + 1; iedge2 < h->edges.end(); iedge2++) {
					pE = *iedge2;
					if(pE->getDestination(NULL)->inTree == true && pE->getSource(NULL)->inTree == true) {
						if(!h->edges.empty()) {
							h->edges.erase(iedge2);
						}
					}
				}
					//Update Destination
				h = hubs[pEdge->getDestination(NULL)->data - 1];
				for(iedge2 = h->edges.begin() + 1; iedge2 < h->edges.end(); iedge2++) {
					pE = *iedge2;
					if(pE->getDestination(NULL)->inTree == true && pE->getSource(NULL)->inTree == true) {
						if(!h->edges.empty()) {
							h->edges.erase(iedge2);
						}
					}
				}
			}
			//  Update Heap
			heap->updateHeap();
		} 
		else {
			//	C is empty
			for (unsigned int j = 0; j < 5*g->getCount(); j++) {
				if (v.empty()) {
					break;
				}
				c.push_back(v.back());
				v.pop_back();
			}
			sort(c.begin(), c.end(), des_cmp_cost);
		}
	}
//  Now that we have all the hubs we need to connect them.
	sort(possConn.begin(), possConn.end(), asc_cmp_plevel);
	while(treeCount != MAX_TREE_SIZE && !possConn.empty()) {
		pEdge = possConn.back();
		if(pEdge->getDestination(NULL)->isConn != true && pEdge->getSource(NULL)->isConn != true) {
			pEdge->getDestination(NULL)->isConn = true;
			pEdge->getSource(NULL)->isConn = true;
			tree.push_back(pEdge);
			treeCount++;
		}
		possConn.pop_back();
	}
	//  Return the degree constrained minimum spanning tree
	return tree;
}

void abdbmst::move(Ant *a) {
	Vertex* vertWalkPtr;
	vertWalkPtr = a->location;
	Edge* edgeWalkPtr;
	int numMoves = 0;
	vector<Edge*>::iterator e;
	double sum = 0.0;
	vector<Range> edges;
	double value;
	Range* current;
	vector<int> v = *a->visited;
	//	Determine Ranges for each edge
	for ( e = vertWalkPtr->edges.begin() ; e < vertWalkPtr->edges.end(); e++ ) {
		edgeWalkPtr = *e;
		Range r;
		r.assocEdge = edgeWalkPtr;
		r.low = sum;
		sum += edgeWalkPtr->pLevel + g->getVerticeWeight(edgeWalkPtr->getDestination(vertWalkPtr)); // changed to include destination weight
		r.high = sum;
		edges.push_back(r);
	}
	while (numMoves < 5) {
		//	Select an edge at random and proportional to its pheremone level
		value = fmod(rand(),(sum+1));
		for (unsigned int i = 0; i < edges.size(); i++) {
			current = &edges[i];
			if (value >= current->low && value < current->high) {
				//	We will use this edge
				edgeWalkPtr = current->assocEdge;
				break;
			}
		}
		//	We have a randomly selected edge, if that edges hasnt already been visited by this ant
		if (v[edgeWalkPtr->getDestination(vertWalkPtr)->data] == 0) {
			edgeWalkPtr->pUpdatesNeeded++;
			a->location = edgeWalkPtr->getDestination(vertWalkPtr);
			v[edgeWalkPtr->getDestination(vertWalkPtr)->data] = 1; 
			break;
		} else {
			numMoves++;
		}
	}
}

#endif