//  Author: Christopher Lee Jackson
//  

#include <limits>
#include <iostream>
#include <iomanip>
#include <vector>
#include "exceptions.h"

using namespace std;

class BinaryHeap {
public:
	explicit BinaryHeap (int capacity = 100);
	explicit BinaryHeap (const vector<int> &items );
	
	bool isEmpty() const;
	vector<int> deleteMax();
    void updateHeap();
private:
	int currentSize; // Number of elements in heap
	vector<int> array; // the heap array
	void buildHeap();
	void percolateDown(int hole);
};

bool BinaryHeap::isEmpty() const{
    if(currentSize == 0)
        return true;
    return false;
}

vector<int> BinaryHeap::deleteMax() {
    vector<int>* temp;
    if( isEmpty())
        throw UnderflowException();
    temp = array[1];
    array[1] = array[ currentSize-- ];
    percolateDown(1);
    return temp;
}

void BinaryHeap::percolateDown(int hole) {
    int child;
    int tmp = array[hole];
    for( ; hole * 2 <= currentSize; hole = child) {
        child = hole * 2;
        if( child != currentSize && array[child+1]> array[child])
            child++;
        if( array[child] > tmp)
            array[hole] = array[child];
        else
            break;
    }
    array[hole] = tmp;
}

BinaryHeap::BinaryHeap(const vector<int> & items ) : array(items.size() + 10), currentSize( items.size()) {
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

