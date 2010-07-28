//  Author: Christopher Lee Jackson
//  

#ifndef __BINARYHEAP_H_INCLUDED__  
#define __BINARYHEAP_H_INCLUDED__  

#include <limits>
#include <iostream>
#include <iomanip>
#include <vector>
#include "Graph.h"
#include "exceptions.h"

using namespace std;

class BinaryHeap {
public:
	//explicit BinaryHeap ();
	explicit BinaryHeap (const vector<Hub*> &items );
    BinaryHeap(int x = 100);
	~BinaryHeap(){array.~vector();}
	bool isEmpty() const;
    void insert(Hub* x);
	Hub* deleteMax();
    void updateHeap();
	int topSize();
	bool find(Vertex* x);
private:
	vector<Hub*> array; // the heap array
	unsigned int currentSize; // Number of elements in heap
	void buildHeap();
	void percolateDown(int hole);
};

BinaryHeap::BinaryHeap(int x ) : array(x + 10), currentSize( 0) {
}

bool BinaryHeap::find(Vertex* x){
	vector<Hub*>::iterator h;
	Hub* pHub;
	for (h = array.begin(); h < array.end(); h++){
		pHub = *h;
		if(pHub->vertId == x->data)
			return true;
	}
	return false;

}
void BinaryHeap::insert(Hub* x) {
    if(currentSize == array.size() - 1)
        array.resize(array.size() * 2);
    //  Percolate up
    int hole = ++currentSize;
    for(; hole > 1 && x->edges.size() > array[hole /2]->edges.size(); hole /= 2) {
        array[hole] = array[hole/2];
    }
    array[hole] = x;
}

bool BinaryHeap::isEmpty() const{
    if(currentSize == 0)
        return true;
    return false;
}

int BinaryHeap::topSize() {
	return array[1]->edges.size();
}

Hub* BinaryHeap::deleteMax() {
    Hub* temp = new Hub();
    if( isEmpty())
        throw UnderflowException();
    temp = array[1];
    array[1] = array[ currentSize-- ];
    percolateDown(1);
    return temp;
}

void BinaryHeap::percolateDown(int hole) {
    unsigned int child;
    Hub* tmp = array[hole];
    for( ; hole * (unsigned int) 2 <= currentSize; hole = child) {
        child = hole * 2;
        if( child != currentSize && array[child+1]->edges.size() > array[child]->edges.size())
            child++;
        if( array[child]->edges.size() > tmp->edges.size())
            array[hole] = array[child];
        else
            break;
    }
    array[hole] = tmp;
}

BinaryHeap::BinaryHeap(const vector<Hub*> & items ) : array(items.size() + 10), currentSize( items.size()) {
    for( unsigned int i = 0; i < items.size(); i++)
        array[i+1] = items[i];
    buildHeap();
}

void BinaryHeap::updateHeap() {
    buildHeap();
}

void BinaryHeap::buildHeap() {
    for(int i = currentSize / 2; i > 0; i--)
        percolateDown(i);
}

#endif
