//  Author: Christopher Lee Jackson
//  Course: CMPMSC 463 - Spring 10


#ifndef __GRAPH_H_INCLUDED__  
#define __GRAPH_H_INCLUDED__  


#include <limits>
#include <iostream>
#include <vector>
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
    ~Hub(){}
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
	bool inTree;
    bool isConn;
	
	double x_coord, y_coord;
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
    double weight;
	double pLevel;
	void setSource(Vertex* s);
	void setDestination(Vertex* d);
	Vertex* getSource(Vertex* loc);
	Vertex* getDestination(Vertex* loc);
};

void Edge::setSource(Vertex* s) {
	source = s;
}

void Edge::setDestination(Vertex* d) {
	destination = d;
}

Vertex* Edge::getSource(Vertex* loc) {
	if(!loc)
		return source;
	else if( loc->data == source->data)
		return source;
	else
		return destination;
}


Vertex* Edge::getDestination(Vertex* loc) {
	if(!loc)
		return destination;
	else if( loc->data == source->data)
		return destination;
	else
		return source;
}	//	END EDGE


/*
 *
 *	GRAPH CLASS
 *
 */

class Graph{
private:
    unsigned int count;
    Vertex *first;
	
public:
    Graph();
    ~Graph();
    int insertVertex(int dataIn);
    int deleteVertex(int dltKey);
    int insertEdge (int fromKey, int toKey, double weight);
    double insertEdge(int fromKey, int toKey);
    int insertVertex(int dataIn, double x, double y);
    bool emptyGraph();
    unsigned int getCount();
	void print();
    void search(Vertex *vertPtr);
	Vertex* getFirst();
    double getVerticeWeight(Vertex *vertPtr);
};

Graph::Graph() {
    count = 0;
    first = NULL;
}

Graph::~Graph() {
    
}

double Graph::getVerticeWeight(Vertex *vertPtr) {
    double sum  = 0.0;
    Edge* edgeWalkPtr;
    vector<Edge*>::iterator e;
    for ( e = vertPtr->edges.begin(); e < vertPtr->edges.end(); e++ ) {
		edgeWalkPtr = *e;
        sum += edgeWalkPtr->pLevel;
    }
    return sum;
}

Vertex* Graph::getFirst() {
	return first;
}

bool Graph::emptyGraph() {
    return (count == 0);
}


/*
 Insert data into the graph.
 */
int Graph::insertVertex(int dataIn) {
	Vertex *newPtr;
    Vertex *locPtr;
    Vertex *predPtr;
    newPtr = new Vertex;
    if(newPtr) {
        newPtr->pNextVert = NULL;
        newPtr->data = dataIn;
        newPtr->degree = 0;
        newPtr->visited = 0;
        newPtr->inTree = false;
        newPtr->isConn = false;
       // newPtr->edges = new vector<Edge*>;
        count++;
    } else {
        //  Memory overflow
        return -1;
    }
    locPtr = first;
    if(!locPtr) {
        //  Inserting into empty graph
        first = newPtr;
    } else {
        predPtr = NULL;
        while(locPtr && dataIn > (locPtr->data)) {
            predPtr = locPtr;
            locPtr = locPtr->pNextVert;
        }
        if (!predPtr) {
            //  insert before fist vertex
            first = newPtr;
        } else {
            predPtr->pNextVert = newPtr;
        }
        newPtr->pNextVert = locPtr;
    }
    return 1;
}

int Graph::insertVertex(int dataIn, double x, double y) {
	Vertex *newPtr;
    Vertex *locPtr;
    Vertex *predPtr;
    newPtr = new Vertex;
    if(newPtr) {
        newPtr->pNextVert = NULL;
        newPtr->data = dataIn;
        newPtr->degree = 0;
        newPtr->visited = 0;
        newPtr->inTree = false;
        newPtr->isConn = false;
        newPtr->x_coord = x;
        newPtr->y_coord = y;
       // newPtr->edges = new vector<Edge*>;
        count++;
    } else {
        //  Memory overflow
        return -1;
    }
    locPtr = first;
    if(!locPtr) {
        //  Inserting into empty graph
        first = newPtr;
    } else {
        predPtr = NULL;
        while(locPtr && dataIn > (locPtr->data)) {
            predPtr = locPtr;
            locPtr = locPtr->pNextVert;
        }
        if (!predPtr) {
            //  insert before fist vertex
            first = newPtr;
        } else {
            predPtr->pNextVert = newPtr;
        }
        newPtr->pNextVert = locPtr;
    }
    return 1;
}


/*
 Delete an existing vertex only if its degree is 0.
 */
int Graph::deleteVertex(int dltKey) {
    Vertex *predPtr;
    Vertex *walkPtr;
    if(!first) {
        return -2;
    }
    predPtr = NULL;
    walkPtr = first;
    while(walkPtr && dltKey > (walkPtr->data)) {
        predPtr = walkPtr;
        walkPtr = walkPtr->pNextVert;
    }
    if (!walkPtr || dltKey != (walkPtr->data)) {
        return -2;
    }
    //  Found vertex, test its degree
    if(walkPtr->degree > 0) {
        //  Cant delte
        return -1;
    }
    //  Can delete
    if(!predPtr) {
        first = walkPtr->pNextVert;
    } else {
        predPtr->pNextVert = walkPtr->pNextVert;
    }
    count--;
    delete walkPtr;
    return 1;
}

/*
 Insert an edge between two verticies.
 */
int Graph::insertEdge(int fromKey, int toKey, double weight) {
    Edge *newPtr;

    Vertex *vertFromPtr;
    Vertex *vertToPtr;
    
    newPtr = new Edge;
    newPtr->weight = weight;
    if(!newPtr) {
        return (-1);
    }
    //  Find source vertex
    vertFromPtr = first;
    while(vertFromPtr && fromKey > (vertFromPtr->data)) {
        vertFromPtr = vertFromPtr->pNextVert;
    }
    if(!vertFromPtr || fromKey != (vertFromPtr->data)) {
        return (-2);
    }
    //  Find destination vertex
    vertToPtr = first;
    while(vertToPtr && toKey > (vertToPtr->data)) {
        vertToPtr = vertToPtr->pNextVert;
    }
    if(!vertToPtr || toKey != (vertToPtr->data)) {
        return (-3);
    }
    //  Found verticies. Make edge.
    ++vertFromPtr->degree;
    ++vertToPtr->degree;
    newPtr->setDestination(vertToPtr);
	newPtr->setSource(vertFromPtr);
	//	Add edges to each adjacency list
	vertToPtr->edges.push_back(newPtr);
	vertFromPtr->edges.push_back(newPtr);
	return 1;
}

double Graph::insertEdge(int fromKey, int toKey) {
    Edge *newPtr;
	
	double weight = 0;
	
    Vertex *vertFromPtr;
    Vertex *vertToPtr;
    
    newPtr = new Edge;

    if(!newPtr) {
        return (-1);
    }
    //  Find source vertex
    vertFromPtr = first;
    while(vertFromPtr && fromKey > (vertFromPtr->data)) {
        vertFromPtr = vertFromPtr->pNextVert;
    }
    if(!vertFromPtr || fromKey != (vertFromPtr->data)) {
        return (-2);
    }
    //  Find destination vertex
    vertToPtr = first;
    while(vertToPtr && toKey > (vertToPtr->data)) {
        vertToPtr = vertToPtr->pNextVert;
    }
    if(!vertToPtr || toKey != (vertToPtr->data)) {
        return (-3);
    }
    //  Found verticies. Make edge.
    weight = sqrt((((vertFromPtr->x_coord - vertToPtr->x_coord) * (vertFromPtr->x_coord - vertToPtr->x_coord)) 
    	+ ((vertFromPtr->y_coord - vertToPtr->y_coord) * (vertFromPtr->y_coord - vertToPtr->y_coord))));
    newPtr->weight = weight;
    ++vertFromPtr->degree;
    ++vertToPtr->degree;
    newPtr->setDestination(vertToPtr);
	newPtr->setSource(vertFromPtr);
	//	Add edges to each adjacency list
	vertToPtr->edges.push_back(newPtr);
	vertFromPtr->edges.push_back(newPtr);
	return weight;
}


unsigned int Graph::getCount() {
    return count;
}

void Graph::print() {
    Vertex *vertWalkPtr;
    //  Set all nodes to unvisited
    vertWalkPtr = first;
    while(vertWalkPtr) {
        vertWalkPtr->visited = 0;
        vertWalkPtr = vertWalkPtr->pNextVert;
    }
    //  Process each vertex in list
    vertWalkPtr = first;
    while(vertWalkPtr) {
        if(vertWalkPtr->visited == 0) {
            search(vertWalkPtr);
        }
        vertWalkPtr = vertWalkPtr->pNextVert;
    }
}


void Graph::search(Vertex *vertPtr) {
    //  Set vertex to processed
    vertPtr->visited =1;
	cout << "Vertex: " << vertPtr->data << ", has edges to: " << endl;
	vector<Edge*>::iterator e;
	for ( e = vertPtr->edges.begin() ; e < vertPtr->edges.end(); e++ ) {
		Edge* c = *e;
		cout << "Vertex: " << c->getDestination(vertPtr)->data << ", with cost: " << c->weight << ", with ph: " << c->pLevel << endl;
	}
	cout << "//end of this vertex" << endl << endl;
    //  Check each Edge for the vertex
	for ( e = vertPtr->edges.begin() ; e < vertPtr->edges.end(); e++ ) {
		Edge* d = *e;
		if (d->getDestination(vertPtr)->visited == 0) {
			search(d->getDestination(vertPtr));
		}
	}
}	//	END GRAPH

#endif
