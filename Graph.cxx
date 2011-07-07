#include "Graph.h"
#include <limits>
#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <cstdlib>

/************************************
***                               ***
***       Edge Definitions        ***
***                               ***
************************************/

void Edge::setSource(Vertex* s) {
    a = s;
}

void Edge::setDestination(Vertex* d) {
    b = d;
}

Vertex* Edge::getSource(Vertex* loc) {
    if(!loc)
        return a;
    else if( loc->data == a->data)
        return b;
    else
        return a;
}

Vertex* Edge::getOtherSide(Vertex* loc) {
    if(loc->data == a->data)
        return b;
    else
        return a;
}

Vertex* Edge::getDestination(Vertex* loc) {
    if(!loc)
        return b;
    else if( loc->data == a->data)
        return b;
    else
        return a;
}   //  END EDGE

void Vertex::updateVerticeWeight() {
    sum  = 0.0;
    Edge* edgeWalkPtr;
    vector<Edge*>::iterator e;
    for ( e = edges.begin(); e < edges.end(); e++ ) {
        edgeWalkPtr = *e;
        sum += edgeWalkPtr->pLevel;
    }
}

/************************************
***                               ***
***       Graph Definitions       ***
***                               ***
************************************/


Graph::Graph() {
    numNodes = 0;
    root = 0;
    oddRoot = 0;
    first = NULL;
    vDepths = new vector<Vertex*>*[50];
    for( int i =0 ; i < 50 ; i++ )
        vDepths[i] = new vector<Vertex*>;
}

Graph::~Graph() {
    for( int i =0 ; i < 50 ; i++ )
        delete vDepths[i];
    delete []vDepths;
}

Vertex* Graph::getFirst() {
    return first;
}

Vertex* Graph::getRand() { 
    int x = rg.IRandom(0,((numNodes - 1)));
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
        newPtr->sum = 0.0;
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
    nodes[dataIn] = newPtr;
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
        newPtr->sum = 0.0;
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

    nodes[dataIn] = newPtr;
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
    newPtr->b=vertToPtr;
    newPtr->a=vertFromPtr;
    //  Add edges to each adjacency list
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
    newPtr->b=vertToPtr;
    newPtr->a=vertFromPtr;
    //  Add edges to each adjacency list
    vertToPtr->edges.push_back(newPtr);
    vertFromPtr->edges.push_back(newPtr);
    return weight;
}

/*
 * Remove a given edge from the graph
 */
void Graph::removeEdge(int a, int b){
    vector<Edge*>::iterator e;
    Edge* eWalkPtr;
    vector<Edge*>::iterator end = nodes[a]->edges.end();
    for ( e = nodes[a]->edges.begin() ; e < end; e++ ) {
        eWalkPtr = *e;
        if(eWalkPtr->a->data == a && eWalkPtr->b->data == b){
            //removeEdge from vector
            nodes[a]->edges.erase(e);
            break;
        }
        else if(eWalkPtr->a->data == b && eWalkPtr->b->data == a){
            //removeEdge from vector
            nodes[a]->edges.erase(e);
            break;

        }
    }
    end = nodes[b]->edges.end();
    for ( e = nodes[b]->edges.begin() ; e < end; e++ ) {
        eWalkPtr = *e;
        if(eWalkPtr->a->data == a && eWalkPtr->b->data == b){
            //removeEdge from vector
            nodes[b]->edges.erase(e);
            break;
        }
        else if(eWalkPtr->a->data == b && eWalkPtr->b->data == a){
            //removeEdge from vector
            nodes[b]->edges.erase(e);
            break;
        }
    }
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
    int unsigned x = 0;
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
    x++;
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
                x++;
                eWalkPtr->getOtherSide(vertWalkPtr)->visited = true;
            }
            //cout << "blah" << endl;
        }
    }
    if(x != numNodes)
        return -1;
    return i;
}

Vertex* Graph::BFS_2(Vertex* pVert) {
    Edge *eWalkPtr;
    vector<Edge*>::iterator e;
    queue<Vertex*> q;
    Vertex *vertWalkPtr;
    Vertex *otherSide;
    vertWalkPtr = first;
    while(vertWalkPtr) {
        vertWalkPtr->visited = false;
        vertWalkPtr = vertWalkPtr->pNextVert;
    }
    q.push(pVert);
    pVert->visited = true;
    pVert->depth = 0;
    vDepths[pVert->depth]->push_back(pVert);
    while(!q.empty()) {
        vertWalkPtr = q.front();
        q.pop();
       
        for ( e = vertWalkPtr->edges.begin() ; e < vertWalkPtr->edges.end(); e++ ) {
            eWalkPtr = *e;
            otherSide = eWalkPtr->getOtherSide(vertWalkPtr);
            if(otherSide->visited == false) {
                q.push(otherSide);
                otherSide->visited = true;
                if(otherSide->data == oddRoot) {
                    otherSide->depth = 0;

                }
                else {
                    otherSide->depth = vertWalkPtr->depth + 1;                    
                }
                vDepths[otherSide->depth]->push_back(otherSide);
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
}   //  END GRAPH


int Graph::testDiameter() {
    int max = 0;
    Vertex* pVert;
    pVert = this->BFS_2(nodes[root]);
    max = this->BFS(pVert);
    return max;
}

bool Graph::isConnected() {
	//	This function will test if the graph is connected. It will only work if testDiameter has just been 
	//	called as it relies on the setting of the visited field in each vertex. This was done to save processing
	//	time. 
	bool connected = true;
	Vertex* vertWalkPtr = first;
	while(vertWalkPtr) {
		if(vertWalkPtr->visited == false)
			connected = false;
		vertWalkPtr = vertWalkPtr->pNextVert;
	}
	return connected;
}
