#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include "Graph.h"
using namespace std;

void processFile(Graph *, char* fileName);

double maxCost = 0;
double minCost = std::numeric_limits<double>::infinity();

int main(int argc, char *argv[])
{
    //  Process input from command line
    if (argc != 2) {
        cerr << "Wrong input.\n";
    }
    char* fileName = new char[50];
    strcpy(fileName,argv[1]);
    //  Process input file and get resulting graph
	Graph *g = new Graph();
    processFile(g, fileName);
    g->print();
    return 0;
}


/*void processFile(Graph *g, char* fileName) {
    double x,y,cost;
    //  Open file for reading  
    ifstream inFile;
    inFile.open(fileName);
    assert(inFile.is_open());
    int eCount, vCount;
    //  Create each vertex after getting vertex count
    inFile >> vCount;
    for(int i = 1; i <= vCount; i++) {
    	inFile >> x >> y;
        g->insertVertex(i, x, y);
    }
    //  Create each edge after processing edge count
    eCount = vCount*(vCount-1)/2;
    for(int v1 = 1; v1<= vCount; v1++) {
    	for(int j = 1; v1 + j<= vCount; j++){
        	cost = g->insertEdge(v1, v1 + j);
			if (cost > maxCost)
				maxCost = cost;
			if (cost < minCost)
				minCost = cost;
    	}
	}
}*/

void processFile(Graph *g, char* fileName) {
    double x,y,cost;
    //  Open file for reading  
    ifstream inFile;
    inFile.open(fileName);
    assert(inFile.is_open());
    int eCount, vCount;
    //  Create each vertex after getting vertex count
    inFile >> vCount;
    for(int i = 1; i <= vCount; i++) {
        g->insertVertex(i);
    }
    //  Create each edge after processing edge count
    eCount = vCount*(vCount-1)/2;
    for(int v1 = 1; v1<= vCount; v1++) {
    	for(int j = 1; j <= vCount; j++){
        	inFile >> cost;
        	if(j > v1){
        		g->insertEdge(v1, j, cost);
				if (cost > maxCost)
					maxCost = cost;
				if (cost < minCost)
					minCost = cost;
    		}
   		}
   	}
}
