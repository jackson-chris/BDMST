//  Author: Christopher Lee Jackson
//  

#include <limits>
#include <iostream>
#include <iomanip>
#include <vector>
#include "Graph.h"
#include "exceptions.h"

using namespace std;

class BinaryHeap {
public:
	explicit BinaryHeap (int capacity = 100);
	explicit BinaryHeap (const vector<Hub> &items );
	
	bool isEmpty() const;
	vector<Hub> deleteMax();
    void updateHeap();
private:
	int currentSize; // Number of elements in heap
	vector<Hub> array; // the heap array
	void buildHeap();
	void percolateDown(int hole);
};

bool BinaryHeap::isEmpty() const{
    if(currentSize == 0)
        return true;
    return false;
}

vector<Hub> BinaryHeap::deleteMax() {
    vector<Hub> temp;
    if( isEmpty())
        throw UnderflowException();
    temp = array[1];
    array[1] = array[ currentSize-- ];
    percolateDown(1);
    return temp;
}

void BinaryHeap::percolateDown(int hole) {
    int child;
    Hub tmp = array[hole];
    for( ; hole * 2 <= currentSize; hole = child) {
        child = hole * 2;
        if( child != currentSize && array[child+1]->edges->getCount() > array[child]->edges->getCount())
            child++;
        if( array[child]->edges->getCount() > tmp->edges->getCount())
            array[hole] = array[child];
        else
            break;
    }
    array[hole] = tmp;
}

explicit BinaryHeap::BinaryHeap(const vector<Hub> & items ) : array(items.size() + 10), currentSize( items.size()) {
    for( int i = 0; i < items.size(); i++)
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

