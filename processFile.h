#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <vector>
#include "Graph.h"
#include <cmath>
#include <cstring>
#include <stack>
#include <queue>
#include <cstdlib>



class processFile{
	private:
	double maxCost;
	double minCost;
	public:
	processFile(){maxCost = 0; minCost = std::numeric_limits<double>::infinity();}
	void processFileOld(Graph *g, ifstream &inFile);
	void processEFile(Graph *g, ifstream &inFile);
	void processRFile(Graph *g, ifstream &inFile);
	void reset(){maxCost = 0; minCost = std::numeric_limits<double>::infinity();}
	double getMax(){return maxCost;}
	double getMin(){return minCost;}
	
};


void processFile::processFileOld(Graph *g, ifstream &inFile) {
   
    int eCount, vCount;
    int i, j;
    double cost;
    //  Create each vertex after getting vertex count
    inFile >> vCount;
	cout << vCount << endl;
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

void processFile::processEFile(Graph *g, ifstream &inFile) {
    double x,y,cost;
    int eCount, vCount;
    //  Create each vertex after getting vertex count
    inFile >> vCount;
	cout << vCount << endl;
    for(int i = 1; i <= vCount; i++) {
    	inFile >> x >> y;
        g->insertVertex(i, x, y);
    }
    //  Create each edge after processing edge count
    eCount = vCount*(vCount-1)/2;
    for(int v1 = 1; v1 <= vCount; v1++) {
    	for(int j = 1; v1 + j <= vCount; j++) {
        	cost = g->insertEdge(v1, v1 + j);
			if (cost > maxCost) {
				maxCost = cost;
			}
			if (cost < minCost) {
				minCost = cost;
			}
    	}
	}
}

void processFile::processRFile(Graph *g, ifstream &inFile) {
    double cost;
    int eCount, vCount;
    //  Create each vertex after getting vertex count
    inFile >> vCount;
	cout << vCount << endl;
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
