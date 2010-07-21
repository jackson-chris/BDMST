//	Author: Christopher Lee Jackson & Jason Jones
//	Course: CMPSC463
//  Problem: 3-1
//  Description:


#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <vector>
#include "Graph.h"
#include "abdbmst.h"
#include <cmath>
#include <cstring>

using namespace std;

void processEFile(Graph *g, char* fileName, int d);
void processRFile(Graph *g, char* fileName, int d);
void processFileOld(Graph *g, char* fileName, int d);

int main( int argc, char *argv[])
{
//	Process input from command line
	if (argc != 4) {
		cerr << "Usage: ab_dbmst <input.file> <fileType> <diameter_bound>\n";
		cerr << "Where: fileType: r = random, e = estien or euc, o = other\n";
		cerr << "       diameter_bound: an integer i, s.t. 4 <= i < |v| \n";
		return 1;
	}
	char* fileName = new char[50];
	char* fileType = new char[2];
	strcpy(fileName,argv[1]);
	strcpy(fileType,argv[2]);
	int d;
	d = atoi(argv[3]);
//  Process input file and get resulting graph
	Graph *g;
	if(fileType[0] == 'e') {
		cout << "USING e file type" << endl;
		processEFile(g, fileName, d);
	}
	else if (fileType[0] == 'r') {
		cout << "USING r file type" << endl;
		processRFile(g, fileName, d);
	}
	else {
		cout << "USING o file type" << endl;
		processFileOld(g, fileName, d);
	}
	return 0;
}
 
void processFileOld(Graph *g, char* fileName, int d) {
	//  Open file for reading
	g = new Graph();
	abdbmst* a;
	ifstream inFile;
	inFile.open(fileName);
	assert(inFile.is_open());
	int eCount, vCount;
	int i, j;
	double cost, maxCost = 0, minCost = std::numeric_limits<double>::infinity();
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
	//	Now that we have the graph lets process it
	a = new abdbmst(g,d, minCost, maxCost);
	delete a;
}

void processEFile(Graph *g, char* fileName, int d) {
	abdbmst* a;
	double x,y,cost, maxCost, minCost;
	int numInstances = 0;
	vector<Edge*> best;
	//  Open file for reading  
	ifstream inFile;
	inFile.open(fileName);
	assert(inFile.is_open());
	int eCount, vCount;
	//	Get number of instances
	inFile >> numInstances;
	for( int n = 0; n < numInstances; n++) {
		g = new Graph();
		maxCost = 0;
		minCost = std::numeric_limits<double>::infinity();
		//  Create each vertex after getting vertex count
		inFile >> vCount;
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
		//	Now that we have the graph lets process it
		a = new abdbmst(g,d,minCost,maxCost);
	}

}

void processRFile(Graph *g, char* fileName, int d) {
	abdbmst* a;
	int numInstances = 0;
	double cost, minCost, maxCost;
	vector<Edge*> best;
	//  Open file for reading  
	ifstream inFile;
	inFile.open(fileName);
	assert(inFile.is_open());
	int eCount, vCount;
	//	Get number of instances
	inFile >> numInstances;
	for (int n = 0; n < numInstances; n++) {
		g = new Graph();
		maxCost = 0;
		minCost = std::numeric_limits<double>::infinity();
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
		//	Now that we have the graph lets process it
		a = new abdbmst(g,d,minCost,maxCost);
	}

}