#include "Graph.h"
#include <limits>
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <cstdlib>

/************************************
***								  ***
***       Edge Definitions        ***
***								  ***
************************************/

void Edge::setSource(Vertex* s) {
    a = s;
}

void Edge::setDestination(Vertex* d) {
	b = d;
}

/*Vertex* Edge::getSource(Vertex* loc) {
	if(!loc)
		return source;
	else if( loc->data == source->data)
		return source;
	else
		return destination;
}*/

Vertex* Edge::getOtherSide(Vertex* loc) {
    if(loc->data == a->data)
        return a;
    else
        return b;
}

/*Vertex* Edge::getDestination(Vertex* loc) {
	if(!loc)
		return destination;
	else if( loc->data == source->data)
		return destination;
	else
		return source;
}*/	//	END EDGE

/************************************
***								  ***
***       Graph Definitions       ***
***								  ***
************************************/


Graph::Graph() {
    numNodes = 0;
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

Vertex* Graph::getRand() {
	int x = rand() % numNodes;
	Vertex* randPtr = first;
	for(int i = 0; i < x; i++){
		randPtr = randPtr->pNextVert;
	}
	return randPtr;
}

bool Graph::emptyGraph() {
    return (numNodes == 0);
}


/*
 Insert data into the graph.
 */
int Graph::insertVertex(int dataIn, Hub* hub) {
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
        newPtr->depth = -1;
        newPtr->pHub = hub;
       // newPtr->edges = new vector<Edge*>;
        numNodes++;
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
        newPtr->depth = -1;
        newPtr->x_coord = x;
        newPtr->y_coord = y;
       // newPtr->edges = new vector<Edge*>;
        numNodes++;
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
    numNodes--;
    delete walkPtr;
    return 1;
}

/*
 Insert an edge between two verticies.
 */
int Graph::insertEdge(int fromKey, int toKey, double weight, double level) {
    Edge *newPtr;

    Vertex *vertFromPtr;
    Vertex *vertToPtr;
    
    newPtr = new Edge;
    newPtr->weight = weight;
    newPtr->pLevel = level;
    newPtr->usable = true;
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

unsigned int Graph::getNumNodes() {
    return numNodes;
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
            print_search(vertWalkPtr);
        }
        vertWalkPtr = vertWalkPtr->pNextVert;
    }
}

int Graph::BFS(Vertex* pVert) {
    int i = 0;
    Edge *eWalkPtr;
    vector<Edge*>::iterator e;
    Vertex *b = new Vertex();
    b->data = -1;
    b->visited = true;
    queue<Vertex*> q;
    Vertex *vertWalkPtr;
    vertWalkPtr = first;
    while(vertWalkPtr) {
        vertWalkPtr->visited = false;
        vertWalkPtr = vertWalkPtr->pNextVert;
    }
    q.push(pVert);
    pVert->visited = true;
    q.push(b);
    while(!q.empty()) {
        vertWalkPtr = q.front();
        q.pop();
        if(vertWalkPtr->data == -1 && !q.empty()) {
            i++;
            q.push(b);
        }
        //cout << vertWalkPtr->data << ", i: " << i << endl;
        for ( e = vertWalkPtr->edges.begin() ; e < vertWalkPtr->edges.end(); e++ ) {
            eWalkPtr = *e;
            if(eWalkPtr->getOtherSide(vertWalkPtr)->visited == false) {
              //  cout << "pushing "<< eWalkPtr->getOtherSide(vertWalkPtr) << endl;
                q.push(eWalkPtr->getOtherSide(vertWalkPtr));
                eWalkPtr->getOtherSide(vertWalkPtr)->visited = true;
            }
            //cout << "blah" << endl;
        }
    }
    return i;
}

Vertex* Graph::BFS_2(Vertex* pVert) {
    Edge *eWalkPtr;
    vector<Edge*>::iterator e;
    queue<Vertex*> q;
    Vertex *vertWalkPtr;
    vertWalkPtr = first;
    while(vertWalkPtr) {
        vertWalkPtr->visited = false;
        vertWalkPtr = vertWalkPtr->pNextVert;
    }
    q.push(pVert);
    pVert->visited = true;
    while(!q.empty()) {
        vertWalkPtr = q.front();
        q.pop();
       
        for ( e = vertWalkPtr->edges.begin() ; e < vertWalkPtr->edges.end(); e++ ) {
            eWalkPtr = *e;
            if(eWalkPtr->getOtherSide(vertWalkPtr)->visited == false) {
                q.push(eWalkPtr->getOtherSide(vertWalkPtr));
                eWalkPtr->getOtherSide(vertWalkPtr)->visited = true;
            }
        }
    }
    return vertWalkPtr;
}

void Graph::print_search(Vertex *vertPtr) {
    Edge *c, *d;
    //  Set vertex to processed
    vertPtr->visited =1;
	cout << "Vertex: " << vertPtr->data << ", has edges to: " << endl;
	vector<Edge*>::iterator e;
	for ( e = vertPtr->edges.begin() ; e < vertPtr->edges.end(); e++ ) {
        c = *e;
		cout << "Vertex: " << c->getOtherSide(vertPtr)->data << ", with cost: " << c->weight << ", with ph: " << c->pLevel << endl;
	}
	cout << "//end of this vertex" << endl << endl;
    //  Check each Edge for the vertex
	for ( e = vertPtr->edges.begin() ; e < vertPtr->edges.end(); e++ ) {
        d = *e;
		if (d->getOtherSide(vertPtr)->visited == 0) {
			print_search(d->getOtherSide(vertPtr));
		}
	}
}	//	END GRAPH

