//  File: bdmst.cxx
//  Author: Christopher Lee Jackson & Jason Jones
//  Description: This is our implementation of our ant based algorithm to aproximate the BDMST problem.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <vector>
#include "randomc.h"
#include "mersenne.cxx"

//  Variables for proportional selection
int32 seed = time(0), rand_pher;
//int32 seed = 1310077132;
TRandomMersenne rg(seed);

#include "Graph.cxx"
#include "Queue.h"
#include <cmath>
#include <cstring>
#include <stack>
#include <cstdlib>
#include "processFile.h"

using namespace std;

typedef struct {
    int data; // initial vertex ant started on
    int nonMove;
    Vertex *location;
    Queue* vQueue;
}Ant;

typedef struct {
    double low;
    double high;
    Edge *assocEdge;
}Range;

//  Globals
const double P_UPDATE_EVAP = 1.05;
const double P_UPDATE_ENHA = 1.05;
const unsigned int TABU_MODIFIER = 5;
const int MAX_CYCLES = 2500; // change back to 2500
const int ONE_EDGE_OPT_BOUND = 20; // was 500
const int ONE_EDGE_OPT_MAX = 2500;
const int K = 250;

int instance = 0;
double loopCount = 0;
double evap_factor = .65;
double enha_factor = 1.5;
double maxCost = 0;
double minCost = std::numeric_limits<double>::infinity();

/*//  Variables for proportional selection
int32 seed = time(0), rand_pher;
//int32 seed = 1310007585;
TRandomMersenne rg(seed);*/

int cycles = 1;
int totalCycles = 1;

//  Prototypes
vector<Edge*> AB_DBMST(Graph *g, int d);
vector<Edge*> locOpt(Graph *g, int d, vector<Edge*> *best);
bool asc_cmp_plevel(Edge *a, Edge *b);
bool des_cmp_cost(Edge *a, Edge *b);
bool asc_src(Edge *a, Edge *b);
bool asc_hub(Hub* a, Hub* b);
void move(Graph *g, Ant *a);
void updatePheromonesPerEdge(Graph *g);
void updatePheromonesGlobal(Graph *g, vector<Edge*> *best, bool improved);
void printEdge(Edge* e);
void resetItems(Graph* g, processFile p);
void compute(Graph* g, int d, processFile p, int instance);
bool replenish(vector<Edge*> *c, vector<Edge*> *v, const unsigned int & CAN_SIZE);
Edge* remEdge(vector<Edge*> best);
int find(vector<int> UF, int start);
int universalSearch(Graph* g, int start);
void populateVector(Graph* g, vector<Edge*> *v);
void populateVector_v2(Graph* g, vector<Edge*> *v, int level);
void getCandidateSet(vector<Edge*> *v, vector<Edge*> *c, const unsigned int & CAN_SIZE);
vector<Edge*> opt_one_edge_v1(Graph* g, Graph* gOpt, vector<Edge*> *tree, unsigned int treeCount, int d);
vector<Edge*> opt_one_edge_v2(Graph* g, Graph* gOpt, vector<Edge*> *tree, int d);
vector<Edge*> hope(Graph *g, int d);

int main( int argc, char *argv[]) {
    //  Process input from command line
    if (argc != 4) {
        cerr << "Usage: ab_dbmst <input.file> <fileType> <diameter_bound>\n";
        cerr << "Where: fileType: r = random, e = estien or euc, o = other\n";
        cerr << "       diameter_bound: an integer i, s.t. 4 <= i < |v| \n";
        return 1;
    }
    char* fileName = new char[50];
    char* fileType = new char[2];
    int numInst = 0;
    strcpy(fileName,argv[1]);
    strcpy(fileType,argv[2]);
    int d;
    d = atoi(argv[3]);
    processFile p;
    //  Open file for reading
    ifstream inFile;
    inFile.open(fileName);
    assert(inFile.is_open());
    //  Process input file and get resulting graph
    Graph* g;
    cout << "INFO: Parameters: " << endl;
    cout << "INFO: P_UPDATE_EVAP: " << P_UPDATE_EVAP << ", P_UPDATE_ENHA: " << P_UPDATE_ENHA << ", Tabu_modifier: " << TABU_MODIFIER << endl;
    cout << "INFO: max_cycles: " << MAX_CYCLES << ", evap_factor: " << evap_factor << ", enha_factor: " << enha_factor << endl;
    cout << "INFO: Input file: " << fileName << ", Diameter Constraint: " << d << endl << endl;
    cout << "INFO: MT Seed: " << seed << endl;
    if(fileType[0] == 'e') {
        //cout << "USING e file type" << endl;
        inFile >> numInst;
        cout << "INFO: num_inst: " << numInst << endl;
        cout << "INFO: ";
        for(int i = 0; i < numInst; i++) {
            //cout << "Instance num: " << i+1 << endl;
            g = new Graph();
            p.processEFile(g, inFile);
            instance++;
            compute(g, d, p, i);
            resetItems(g, p);
        }
    }
    else if (fileType[0] == 'r') {
        //cout << "USING r file type" << endl;
        inFile >> numInst;
        for(int i = 0; i < numInst; i++) {
            //cout << "Instance num: " << i+1 << endl;
            g = new Graph();
            p.processRFile(g, inFile);
            compute(g, d, p, i);
            resetItems(g, p);
        }
    }
    else {
        //cout << "USING o file type" << endl;
        g = new Graph();
        p.processFileOld(g, inFile);
        compute(g, d, p, 0);
    }
    return 0;
}

void compute(Graph* g, int d, processFile p, int i) {
    maxCost = p.getMax();
    minCost = p.getMin();
    if((unsigned int) d > g->numNodes) {
        cout << "No need to run this diameter test. Running MST will give you solution, since diameter is greater than number of nodes." << endl;
        exit(1);
    }
    cout << "Instance number: " << i << endl;
    vector<Edge*> best = AB_DBMST(g, d);
    //cout << "Size of best: " << best.size() << endl;
    //sort(best.begin(), best.end(), asc_src);
    //cout << "Best Tree num edges: " << best.size() << endl;
    //for_each(best.begin(), best.end(), printEdge);
}

/*
*
*
*   For Each, Sort - Helper Functions
*
*/
bool asc_cmp_plevel(Edge *a, Edge *b) {
    return (a->pLevel < b->pLevel);
}

bool des_cmp_cost(Edge *a, Edge *b) {
    return (a->weight > b->weight);
}

bool asc_src(Edge* a, Edge* b) {
    return (a->a->data < b->a->data);
}

bool asc_hub(Hub* a, Hub* b) {
    return (a->vertId < b->vertId);
}

void printEdge(Edge* e) {
    cout << "RESULT: ";
    cout << e->a->data << " " << e->b->data << " " << e->weight << " " << e->pLevel << endl;
}

void resetItems(Graph* g, processFile p) {
    free(g);
    maxCost = 0;
    minCost = std::numeric_limits<double>::infinity();
    p.reset();
}

vector<Edge*> AB_DBMST(Graph *g, int d) {
    //  Local Declerations
    double bestCost = std::numeric_limits<double>::infinity();
    double treeCost = 0;
    bool newBest = false;
    const int s = 75;
    int bestRoot = -1, bestOddRoot = -1;
    vector<Edge*> best, current;
    Vertex *vertWalkPtr;
    Edge *edgeWalkPtr, *pEdge;
    vector<Ant*> ants;
    vector<Edge*>::iterator e, ed, iedge1;

    Ant *a;
    //  Assign one ant to each vertex
    vertWalkPtr = g->getFirst();
    for (unsigned int i = 0; i < g->getNumNodes(); i++) {
        a = new Ant;
        a->data = i +1;
        a->nonMove = 0;
        a->location = vertWalkPtr;
        a->vQueue = new Queue(TABU_MODIFIER);
        ants.push_back(a);
        //  Initialize pheremone level of each edge, and set pUdatesNeeded to zero
        for ( e = vertWalkPtr->edges.begin() ; e < vertWalkPtr->edges.end(); e++ ) {
            edgeWalkPtr = *e;
            if (edgeWalkPtr->a == vertWalkPtr) {
                edgeWalkPtr->pUpdatesNeeded = 0;
                edgeWalkPtr->pLevel = (maxCost - edgeWalkPtr->weight) + ((maxCost - minCost) / 3);
            }
        }
        //  Done with this vertex's edges; move on to next vertex
        vertWalkPtr->updateVerticeWeight();
        vertWalkPtr = vertWalkPtr->pNextVert;
    }
    while (totalCycles <= 10000 && cycles <= MAX_CYCLES) { 
        //if(totalCycles % 25 == 0) 
            //cerr << "CYCLE " << totalCycles << endl;
        //  Exploration Stage
        for (int step = 1; step <= s; step++) {
            if (step == s/3 || step == (2*s)/3) {
                updatePheromonesPerEdge(g);
            }
            for (unsigned int j = 0; j < g->getNumNodes(); j++) {
                a = ants[j];
                move(g, a);
            }
        }
        //  Do we even need to still do this since we are using a circular queue?
        for(unsigned int w = 0; w < g->getNumNodes(); w++) {
            ants[w]->vQueue->reset(); //  RESET VISITED FOR EACH ANT
        }
        updatePheromonesPerEdge(g);
        //  Tree Construction Stage
        current = hope(g, d);
        //  Get new tree cost
        for ( ed = current.begin(); ed < current.end(); ed++ ) {
            edgeWalkPtr = *ed;
            treeCost+=edgeWalkPtr->weight;
        }
        if (treeCost < bestCost && (current.size() == g->getNumNodes() - 1)) {
            //cerr << "FOUND NEW BEST at cycle: " << totalCycles <<endl;
            best = current;
            bestCost = treeCost;
            newBest=true;
            bestRoot = g->root;
            bestOddRoot = g->oddRoot;
            if (totalCycles != 1)
                cycles = 0;
        } 
        if (cycles % 100 == 0) {
            if(newBest) {
                updatePheromonesGlobal(g, &best, false);
            } 
            else {
                updatePheromonesGlobal(g, &best, true);
            }
            newBest = false;
        }
        if (totalCycles % 500 == 0) {
            evap_factor *= P_UPDATE_EVAP; 
            enha_factor *= P_UPDATE_ENHA; 
        }
        totalCycles++;
        cycles++;
        treeCost = 0;
    }

    // Test if it meets the diameter bound.
    Graph* gTest = new Graph();
    //  add all vertices
    Vertex* pVert = g->getFirst();
    while(pVert) {
        gTest->insertVertex(pVert->data);
        pVert = pVert->pNextVert;
    }
    //  Now add edges to graph.
    for(iedge1 = best.begin(); iedge1 < best.end(); iedge1++) {
        pEdge = *iedge1;
        gTest->insertEdge(pEdge->a->data, pEdge->b->data, pEdge->weight, pEdge->pLevel);
    }
    gTest->root = bestRoot;
    gTest->oddRoot = bestOddRoot;

    cout << "This is the list of edges BEFORE local optimization: " << endl;
    for_each(best.begin(), best.end(), printEdge);
    cout << "RESULT" << instance << ": Cost: " << bestCost << endl;
    cout << "RESULT: Diameter: " << gTest->testDiameter() << endl;
    //best = opt_one_edge_v1(g, gTest, &best, best.size(), d);
    best = opt_one_edge_v2(g, gTest, &best, d);
    cout << "This is the list of edges AFTER local optimization: " << endl;
    for_each(best.begin(), best.end(), printEdge);
    bestCost = 0;
    for ( ed = best.begin(); ed < best.end(); ed++ ) {
        edgeWalkPtr = *ed;
        bestCost+=edgeWalkPtr->weight;
    }
    cout << "RESULT" << instance << ": Cost: " << bestCost << endl;
    cout << "RESULT: Diameter: " << gTest->testDiameter() << endl;
    //  Reset items
    ants.clear();
    cycles = 1;
    totalCycles = 1;
    current.clear();
    treeCost = 0;
    bestCost = 0;
    //  Return best tree
    return best;
}

void updatePheromonesGlobal(Graph *g, vector<Edge*> *best, bool improved) {
    //  Local Variables
    double pMax = 1000*((maxCost - minCost) + (maxCost - minCost) / 3);
    double pMin = (maxCost - minCost)/3;
    Edge *e;
    double XMax = 0.3;
    double XMin = 0.1;
    double rand_evap_factor;
    double IP;
    vector<Edge*>::iterator ex;

    //  For each edge in the best tree update pheromone levels
    for ( ex = best->begin() ; ex < best->end(); ex++ ) {
        e = *ex;
        IP = (maxCost - e->weight) + ((maxCost - minCost) / 3);
        if (improved) {
        //  IMPROVEMENT so Apply Enhancement
            e->pLevel = enha_factor*e->pLevel;
        } else {
        //  NO IMPROVEMENTS so Apply Evaporation
            rand_evap_factor = XMin + rg.BRandom() * (XMax - XMin) / RAND_MAX;
            e->pLevel = rand_evap_factor*e->pLevel;
        }
        //  Check if fell below minCost or went above maxCost
        if (e->pLevel > pMax) {
            e->pLevel = pMax - IP;
        } else if (e->pLevel < pMin) {
            e->pLevel = pMin + IP;
        }
    }
}


void updatePheromonesPerEdge(Graph *g) {
    //  Local Variables
    Vertex *vertWalkPtr = g->getFirst();
    double pMax = 1000*((maxCost - minCost) + (maxCost - minCost) / 3);
    double pMin = (maxCost - minCost)/3;
    double IP;
    vector<Edge*>::iterator ex;
    Edge *edgeWalkPtr;

    while (vertWalkPtr) {
        for ( ex = vertWalkPtr->edges.begin() ; ex < vertWalkPtr->edges.end(); ex++ ) {
            edgeWalkPtr = *ex;
            if (edgeWalkPtr->a == vertWalkPtr) {
                IP = (maxCost - edgeWalkPtr->weight) + ((maxCost - minCost) / 3);
                edgeWalkPtr->pLevel = (1 - evap_factor)*(edgeWalkPtr->pLevel)+(edgeWalkPtr->pUpdatesNeeded * IP);
                if (edgeWalkPtr->pLevel > pMax) {
                    edgeWalkPtr->pLevel = pMax - IP;
                } else if (edgeWalkPtr->pLevel < pMin) {
                    edgeWalkPtr->pLevel = pMin + IP;
                }
                //  Done updating this edge reset multiplier
                edgeWalkPtr->pUpdatesNeeded = 0;
            }
        }
        vertWalkPtr->updateVerticeWeight();
        vertWalkPtr = vertWalkPtr->pNextVert;
    }
}

void getCandidateSet(vector<Edge*> *v, vector<Edge*> *c, unsigned const int & CAN_SIZE) {
    for (unsigned int i = 0; i < CAN_SIZE; i++) {
        if (v->empty()) {
            break;
        }
        c->push_back(v->back());
        v->pop_back();
    }
    return;
}

void populateVector(Graph* g, vector<Edge*> *v) {
    //  Local Variables
    Vertex* vertWalkPtr;
    vector<Edge*>::iterator ie;
    Edge* edgeWalkPtr;
    //  Logic
    vertWalkPtr = g->getFirst();
    while (vertWalkPtr) {
        vertWalkPtr->treeDegree = 0;
        vertWalkPtr->inTree = false;
        vertWalkPtr->isConn = false;
        for ( ie = vertWalkPtr->edges.begin() ; ie < vertWalkPtr->edges.end(); ie++ ) {
            edgeWalkPtr = *ie;
            //  Dont want duplicate edges in listing
            if (edgeWalkPtr->a == vertWalkPtr) {
                edgeWalkPtr->inTree = false;
                v->push_back(edgeWalkPtr);
            }
        }
        vertWalkPtr = vertWalkPtr->pNextVert;
    }
}

void populateVector_v2(Graph* g, vector<Edge*> *v, int level) {
    //  Local Variables
    Vertex* vertWalkPtr;
    vector<Edge*>::iterator ie;
    vector<Vertex*>::iterator vi;
    Edge* edgeWalkPtr;
    //  Logic
    for ( vi = g->vDepths[level]->begin(); vi < g->vDepths[level]->end(); vi++) {
        vertWalkPtr = *vi;
        vertWalkPtr->inTree = false;
        vertWalkPtr->isConn = false;
        for ( ie = vertWalkPtr->edges.begin() ; ie < vertWalkPtr->edges.end(); ie++ ) {
            edgeWalkPtr = *ie;
            //  Dont want duplicate edges in listing
            if (edgeWalkPtr->getOtherSide(vertWalkPtr)->depth >= level) {
                edgeWalkPtr->inTree = false;
                v->push_back(edgeWalkPtr);
            }
        }
    }
}

vector<Edge*> opt_one_edge_v1(Graph* g, Graph* gOpt, vector<Edge*> *tree, unsigned int treeCount, int d) {
    int rMade = 0;
    Edge* edgeWalkPtr = NULL,* edgeTemp;
    vector<Edge*> newTree;
    int noImp = 0, tries = 0;
    vector<Edge*>::iterator e;
    double sum = 0.0;
    bool improved = false;
    Range* ranges[treeCount];
    int value;
    int i, q;
    int bsint = 0;
    Range* current;
    Range* temp;
    vector<Edge*> v;
    populateVector(g, &v);
    int diameter = 0;
    int numEdge = v.size();
    //initialize ranges
    for ( e = tree->begin(), q=0; e < tree->end(); e++, q++) {
        edgeWalkPtr = *e;
        Range* r = new Range();
        r->assocEdge = edgeWalkPtr;
        r->low = sum;
        sum += edgeWalkPtr->weight * 10000;
        r->high = sum;
        ranges[q] = r;
    }
    //  mark edges already in tree
    for ( e = tree->begin(); e < tree->end(); e++) {
        edgeWalkPtr = *e;
        edgeWalkPtr->inTree = true;
    }

    while (noImp < ONE_EDGE_OPT_BOUND ) {//} && tries < ONE_EDGE_OPT_MAX) {
    //  Pick an edge to remove at random favoring edges with low pheremones
    //  First we determine the ranges for each edge
        value = rg.IRandom(0,((int) (sum))); // produce a random number between 0 and highest range
        i = treeCount / 2;
        if (i%2 != 0)
            i++;
        bsint = i;
        while(true) {
        //cout << "oh shit " << i << "treeCount " << treeCount << endl  ;
            current = ranges[i];
            bsint -= bsint/2;
            if(value < current->low){
                i -= bsint;
            }
            else if(value >= current->high){
                i += bsint;
            }
            else{
            //  We will use this edge
            //cout << current->assocEdge->weight << endl;
                edgeWalkPtr = current->assocEdge;
            break;
            }
        }
        //  We now have an edge that we wish to remove.
        //  Remove the edge
        gOpt->removeEdge(edgeWalkPtr->a->data, edgeWalkPtr->b->data);
        //  Try adding new edge if it improves the tree and doesn't violate the diameter constraint keep it.
        for (int j=0; j < K; j++) {
        //  select a random edge, if its weight is less than the edge we just removed use it to try and improve tree.
            value = rg.IRandom(0, numEdge - 1);
            if (v[value]->weight < edgeWalkPtr->weight && v[value]->inTree == false ) {
                gOpt->insertEdge(v[value]->a->data, v[value]->b->data, v[value]->weight, v[value]->pLevel);
                diameter = gOpt->testDiameter();
                //cout << "the diameter after the addition is: " << diameter << endl;
                if (diameter > 0 && diameter <= d && gOpt->isConnected()) {
                    //cout << "IMPROVEMENT! Lets replace the edge." << endl;
                    edgeWalkPtr->inTree = false;
                    tree->erase(tree->begin() + i);
                    tree->push_back(v[value]);
                    v[value]->inTree = true;
                    // update the associated edge for the range that is being changed.
                    ranges[i]->assocEdge = v[value];
                    //  Update the ranges now that we have added a new edge into the tree
                    if(i != 0)
                        sum = ranges[i-1]->high;
                    else
                        sum = 0;
                    for ( e = tree->begin() + i, q = i ; e < tree->end(); e++, q++ ) {
                        edgeTemp = *e;
                        temp = ranges[q];
                        temp->low = sum;
                        sum += edgeTemp->weight * 10000;
                        temp->high = sum;
                    }
                    improved = true;
                    break;
                } 
                else { 
                    gOpt->removeEdge(v[value]->a->data, v[value]->b->data); 
                }
                if (improved) {
                    break;
                }
            }
            //cout << "END FOR\n";
        }
        //cout << "i broke.\n";
        //  Handle Counters
        if (improved) {
            noImp = 0;
            improved = false;
            rMade++;
        }
        else {
            noImp++;
            gOpt->insertEdge(edgeWalkPtr->a->data, edgeWalkPtr->b->data, edgeWalkPtr->weight, edgeWalkPtr->pLevel);
        }
            tries++;
    }
    //cout << "RESULT: Diameter: " << gOpt->testDiameter() << endl;
    //gOpt->print();
    printf("%d edges were exchanged using opt_one_edge_v1.\n", rMade);
    populateVector(gOpt,&newTree);
    return newTree;
}

vector<Edge*> opt_one_edge_v2(Graph* g, Graph* gOpt, vector<Edge*> *tree, int d) {
    Edge* edgeWalkPtr = NULL, *ePtr = NULL;
    vector<Edge*> newTree, possEdges;
    int tries;
    int levelRemove;
    int levelAdd;
    vector<Edge*>::iterator e;
    double sum;
    int value, updates;
    Range* rWalk;
    int q; 
    vector<Edge*> levelEdges;
    Vertex* vertWalkPtr;
    Range** ranges;
    int tabu_size = (int)(g->numNodes*.10);
    Queue* tQueue = new Queue(tabu_size);


    for(int l = 0; l < ONE_EDGE_OPT_BOUND; l++) {
        levelRemove = rg.IRandom(1,((int) (gOpt->height - 1)));
        levelAdd = levelRemove;
        sum = 0.0;
        value = 0;
        updates = 0;
        tries = 0;
        populateVector_v2(gOpt, &levelEdges, levelRemove);
        //initialize ranges
        ranges = new Range*[levelEdges.size()];
        for(unsigned int k = 0; k < levelEdges.size(); k++) {
            ranges[k] = new Range();
        }
        //cout << "Root" << gOpt->root << endl;
        //cout << "Odd Root " << gOpt->oddRoot << endl;

        if(levelEdges.size() == 0){
            cout << "woops." << endl;
            //populateVector(gOpt,&newTree);
            delete []ranges;
            levelEdges.clear();
            possEdges.clear();
            continue;
            //return newTree;
        }
        for ( e = levelEdges.begin(), q=0; e < levelEdges.end(); e++, q++) {
            edgeWalkPtr = *e;
            Range* r = new Range();
            r->assocEdge = edgeWalkPtr;
            r->low = sum;
            sum += edgeWalkPtr->weight * 10000;
            r->high = sum;
            ranges[q] = r;
        }
        //  mark edges already in tree
        for ( e = tree->begin(); e < tree->end(); e++) {
            edgeWalkPtr = *e;
            edgeWalkPtr->inTree = true;
            edgeWalkPtr->usable = true;
        }
        //cout << levelEdges.size() << endl;
        while (tries < ONE_EDGE_OPT_MAX && updates < 30) {
        //  Pick an edge to remove at random favoring edges with low pheremones
        //  First we determine the ranges for each edge
            edgeWalkPtr = NULL;
            possEdges.clear();
            value = rg.IRandom(0,((int) (sum))); // produce a random number between 0 and highest range
            for(unsigned int i = 0; i < levelEdges.size(); i++) {
                rWalk = ranges[i];
                if(rWalk->low <= value && rWalk->high > value && rWalk->assocEdge->usable && !tQueue->exists(rWalk->assocEdge->a->data) && !tQueue->exists(rWalk->assocEdge->b->data)) {
                    edgeWalkPtr=rWalk->assocEdge;
                    edgeWalkPtr->usable = false;
                }
            }
            if(!edgeWalkPtr)
                break; 
         //  We now have an edge that we wish to remove.
         //  Remove the edge
        
        gOpt->removeEdge(edgeWalkPtr->a->data, edgeWalkPtr->b->data);
        //  update tabu list
        tQueue->push(edgeWalkPtr->a->data);
        tQueue->push(edgeWalkPtr->b->data);
        // find out what vertice we have just cut from.
        vertWalkPtr = edgeWalkPtr->a->depth > edgeWalkPtr->b->depth ? g->nodes[edgeWalkPtr->a->data] : g->nodes[edgeWalkPtr->b->data];
        // Noww get all possible edges for that vertex
        for( e = vertWalkPtr->edges.begin(); e < vertWalkPtr->edges.end(); e++) {
            ePtr = *e;
            if(gOpt->nodes[ePtr->getOtherSide(vertWalkPtr)->data]->depth <= levelAdd)
                possEdges.push_back(ePtr);
        }
        sort(possEdges.begin(), possEdges.end(), des_cmp_cost);
        //cout << endl;
        //for_each(possEdges.begin(), possEdges.end(), printEdge);
        //cout << endl;
        ePtr = possEdges.back();
        //for(unsigned int i = 0; i < possEdges.size(); i++)
            //cout << "possEdges: " << possEdges[i]->a->data << ", " << possEdges[i]->b->data << " " << possEdges[i]->inTree << endl;
        while(ePtr->inTree && !possEdges.empty()){
            possEdges.pop_back();
            ePtr = possEdges.back();
        }
        //cout << "Old Edge" << edgeWalkPtr->a->data << ", " << edgeWalkPtr->b->data << "\t" << edgeWalkPtr->weight << endl;
        //cout << "New Edge" << ePtr->a->data << ", " << ePtr->b->data << "\t" << ePtr->weight << endl;
        //cout << "Depth of new a " << ePtr->a->depth << "Depth of new b " << e
        if(edgeWalkPtr->weight > ePtr->weight) {
            cout << "we improved.\n";
            gOpt->insertEdge(ePtr->a->data, ePtr->b->data, ePtr->weight, ePtr->pLevel);
            updates++;
            break;
        }
        else {
            gOpt->insertEdge(edgeWalkPtr->a->data, edgeWalkPtr->b->data, edgeWalkPtr->weight, edgeWalkPtr->pLevel);
            tries++;
            //cout << "we failed.\n";
        }
        //cout << "try number: " << tries << endl;
    }
    //  reset items
    for(unsigned int k = 0; k < levelEdges.size(); k++) {
        delete ranges[k];
    }
    delete []ranges;
    levelEdges.clear();
    possEdges.clear();
}
populateVector(gOpt,&newTree);
return newTree;
}

vector<Edge*> hope(Graph *g, int d) {
    Vertex* pVert, *pVert2;
    vector<Edge*> v, c, inTree;
    Edge* pEdge;
    vector<Edge*>::reverse_iterator iedge1;
    bool done = false, didReplenish = true;
    //  Put all edges into a vector
    populateVector(g, &v);
    //  Sort edges in ascending order based upon pheromone level
    sort(v.begin(), v.end(), asc_cmp_plevel);
    //  Select 5n edges from the end of v( the highest pheromones edges) and put them into c.
    getCandidateSet(&v,&c,g->numNodes);
    //  Sort edges in descending order based upon cost
    sort(c.begin(), c.end(), des_cmp_cost);
    
    
    pEdge = c.back();
    c.pop_back();
    if ( d%2 == 0 ) {
        g->root = pEdge->a->data;
        pEdge->a->depth=0;
        g->oddRoot = pEdge->b->data;
        pEdge->b->depth=0;
    } else {
        g->root = pEdge->a->data;
        pEdge->a->depth=0;
        pEdge->b->depth=1;
    }
    pEdge->inTree = true;
    pEdge->a->inTree = true;
    pEdge->b->inTree = true;
    inTree.push_back(pEdge);
    while(!done) {
        //cout << "we be looping\n";
        if ( c.empty() ) {
            didReplenish = replenish(&c, &v, g->numNodes);
            //cout << "Replenished" << endl;
        }
        if (!didReplenish)
            break;
        //for(iedge1 = c.end(); iedge1 > c.begin(); iedge1--) {
        for(iedge1 = c.rbegin(); iedge1 < c.rend(); iedge1++) {
            //cout << "yay" << endl;
            pEdge = *iedge1;
            if(pEdge->a->inTree ^ pEdge->b->inTree) {
                pEdge->a->inTree == true ? pVert = pEdge->a : pVert = pEdge->b;
                if (pVert->depth < (d/2) - 1) {
                    //  Add this edge into the tree.
                    pEdge->inTree = true;
                    pVert2 = pEdge->getOtherSide(pVert);
                    pVert2->inTree = true;
                    pVert2->depth = pVert->depth + 1;
                    inTree.push_back(pEdge);
                    //cout << "we added\n";
                    break;
                }
                c.pop_back();
                if(inTree.size() == g->numNodes - 1)
                    done=true;
            }
            if ( pEdge == c.front() ) {
                c.clear();
                //cout << "We cleared." << endl;
            }
        }
    }
    return inTree;
}


int find(vector<int> UF, int start){
    while(UF[start] != start)
        start = UF[start];

    return start;
}

bool replenish(vector<Edge*> *c, vector<Edge*> *v, const unsigned int & CAN_SIZE) {
    if(v->empty()) {
        return false;
    }
    for (unsigned int j = 0; j < CAN_SIZE; j++) {
        if (v->empty()) {
            break;
        }
        c->push_back(v->back());
        v->pop_back();
    }
    sort(c->begin(), c->end(), des_cmp_cost);
    return true;
}

void move(Graph *g, Ant *a) {
    Vertex* vertWalkPtr;
    Vertex* vDest;
    vertWalkPtr = a->location;
    Edge* edgeWalkPtr = NULL;
    int numMoves = 0;
    int size = 0, initialI = 0;
    bool alreadyVisited;
    vector<Edge*>::iterator e;
    double sum = 0.0;
    vector<Range> edges;
    int value;
    int i = 0, bsint = 0;
    Range* current;
    //  Determine Ranges for each edge
    for ( e = vertWalkPtr->edges.begin() ; e < vertWalkPtr->edges.end(); e++ ) {
        edgeWalkPtr = *e;
        Range r;
        r.assocEdge = edgeWalkPtr;
        //cout << "Edge " << ++i << endl;
        r.low = sum;
        //cout << "Low sum " << sum << endl;
        sum += edgeWalkPtr->pLevel;// + edgeWalkPtr->getOtherSide(vertWalkPtr)->sum; 
        r.high = sum;
        //cout << "High sum " << sum << endl << endl << endl;
        edges.push_back(r);
    }
    size = edges.size();
    initialI = size / 2 - 1;
    while (numMoves < 5) {
        //  Select an edge at random and proportional to its pheremone level
        value = rg.IRandom(0,((int) (sum)));
        i = initialI;
        if(i%2 != 0)
            i++;
        bsint = i;
        while(true){
            current = &edges[i];
            bsint -= bsint/2;
            //cout << value << ", " << current->low << ", " << current->high << endl;
            //cout << "i = " << i << " bsint = " << bsint << endl;
            if(value < current->low){
                i -= bsint;
            }
            else if(value >= current->high){
                i += bsint;
            }
            else{
            //  We will use this edge
                //cout << current->assocEdge->weight << endl;
                edgeWalkPtr = current->assocEdge;
                break;
            } } 
        //  Check to see if the ant is stuck
        if (a->nonMove > 4) {
            a->vQueue->reset();
        }
        //  We have a randomly selected edge, if that edges hasnt already been visited by this ant traverse the edge
        vDest = edgeWalkPtr->getOtherSide(vertWalkPtr);
        alreadyVisited = false;
        for(unsigned int j = 0; j <= TABU_MODIFIER; j++) {
            if( a->vQueue->array[j] == vDest->data ) {
                // This ant has already visited this vertex
                alreadyVisited = true;
                break;
            }
        }
        if (!alreadyVisited) {
            edgeWalkPtr->pUpdatesNeeded++;
            a->location = vDest;
            //  the ant has moved so update where required
            a->nonMove = 0;
            a->vQueue->push(vDest->data);
            break;
        } else {
            // Already been visited, so we didn't make a move.
            a->nonMove++;
            numMoves++;
        }
    }
}
