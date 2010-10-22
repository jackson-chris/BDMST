
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
#include <stack>

using namespace std;

double maxCost = 0;
double minCost = std::numeric_limits<double>::infinity();


int testDiameter(Graph* g);
void processFileOld(Graph* g, ifstream &inFile);

int main(int argc, char *argv[]) {
    //  Process input from command line
    if (argc != 2) {
		cerr << "Usage: tree_diam_test <input.file> \n";
		return 1;
	}
    char* fileName = new char[50];
    int numInst = 0;
    strcpy(fileName,argv[1]);
    //  Open file for reading
    ifstream inFile;
    inFile.open(fileName);
    assert(inFile.is_open());
    //  Process input file and get resulting graph
    Graph* g;
	cout << "USING o file type" << endl;
    g = new Graph();
	processFileOld(g, inFile);
    //g->print();
    cout << "diameter of graph is: " << testDiameter(g) << endl;
    //cout << "Vertex #: " << g->getFirst()->pNextVert->data << endl;
    //cout << g->BFS(g->getFirst()) << endl;
    return 0;
}

void processFileOld(Graph *g, ifstream &inFile) {
   
    int eCount, vCount;
    int i, j;
    double cost;
    //  Create each vertex after getting vertex count
    inFile >> vCount;
    for(int i = 1; i <= vCount; i++) {
        g->insertVertex(i);
    }
    //  Create each edge after processing edge count
    eCount = vCount*(vCount-1)/2;
    for(int e = 0; e < eCount; e++) {
        inFile >> i >> j >> cost;
        g->insertEdge(i, j, cost);
		if (cost > maxCost)
			maxCost = cost;
		if (cost < minCost)
			minCost = cost;
    }
}

int testDiameter(Graph* g) {
    int max = 0, temp = 0;
    vector<Edge*>::iterator iEdge;
    //Edge* pEdge;
    Vertex* pVert;
    pVert = g->getFirst();
    while (pVert) {
        temp = g->BFS(pVert);
        if (max < temp)
            max = temp;
        pVert = pVert->pNextVert;
    }
    return max;
}