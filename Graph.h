//  Author: Christopher Lee Jackson
//  Course: CMPMSC 463 - Spring 10


#ifndef __GRAPH_H_INCLUDED__  
#define __GRAPH_H_INCLUDED__  


#include <limits>
#include <iostream>
#include <vector>
#include <queue>
#include <math.h>
#include <cmath>


using namespace std;

class Vertex;
class Edge;
class Hub;

class Hub {
public:
    unsigned int vertId;
    vector<Edge*> edges;
    Vertex* vert;
	bool inTree;
    int lenPath;
    ~Hub(){edges.~vector();}
};

/*
 *
 *	VERTEX CLASS
 *
 */

class Vertex {
public:
    Vertex *pNextVert;
	vector<Edge*> edges;
    int data;
    int degree;
	int treeDegree;
	int visited;
	int depth;
	bool inTree;
    bool isConn;
    Hub* pHub;
	double x_coord, y_coord;
	double sum;
	void updateVerticeWeight();
	~Vertex(){}
};	//	END VERTEX



/*
 *
 *	EDGE CLASS
 *
 */


class Edge {
private:
	Vertex *source;
	Vertex *destination;
public:
	double pUpdatesNeeded;
	bool inTree;
	bool usable;
    double weight;
	double pLevel;
	void setSource(Vertex* s);
	void setDestination(Vertex* d);
	Vertex* getSource(Vertex* loc);
	Vertex* getDestination(Vertex* loc);
    Vertex* getOtherSide(Vertex* loc);
};



/*
 *
 *	GRAPH CLASS
 *
 */

class Graph{
private:
    unsigned int numNodes;
    Vertex *first;
public:
    Graph();
    ~Graph();
    int insertVertex(int dataIn, Hub* hub = NULL);
    int deleteVertex(int dltKey);
    int insertEdge (int fromKey, int toKey, double weight, double pLevel = 0);
    double insertEdge(int fromKey, int toKey);
    int insertVertex(int dataIn, double x, double y);
    bool emptyGraph();
    unsigned int getNumNodes();
	void print();
    void print_search(Vertex *vertPtr);
	Vertex* getFirst();
	Vertex* getRand();
    int BFS(Vertex* pVert);
    Vertex* BFS_2(Vertex* pVert);

};

#endif
