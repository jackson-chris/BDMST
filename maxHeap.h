//  Author: Christopher Lee Jackson
//  Course: CMPMSC 463 - Spring 10

#include <limits>
#include <iostream>
#include <vector>
#include "Graph.h"

using namespace std;

class maxHeap {
private:
	Type heap[] // to do
public:
	void reheapUp(int heap[], int newNode);
	void reheapDown(int heap[], int root, int last);
	bool insertHeap(int& heap[], int& last, int data);
	bool deleteHeap(int heap[], int& last, int& dataOut);
};	//	END maxHeap

void reheapUp(int heap[], int newNode) {
	int parent;
	int hold;
	
	if(newNode) {
		parent = (newNode - 1) / 2;
		if(heap[newNode] > heap[parent]) {
			// child is greater than parent
			hold = heap[parent];
			heap[parent] = heap[newNode];
			heap[newNode] = hold;
			reheapUp(heap,parent);
		}
	}
	return;
}

void reheapDown(int heap[], int root, int last) {
	int hold, leftKey, rightKey, largeChildKey, largeChildIndex;
	if((root * 2 + 1) <= last) {
		// there is at least one child
		leftkey = heap[root*2+1];
		if((root * 2 +2) <= last) 
			rightKey = heap[root * 2 + 2];
		else
			rightKey = -1;
		//Determine which child is larger
		if(leftKey > rightKey) {
			largeChildKey = leftKey;
			largeChildIndex = root * 2 + 1;
		} else {
			largeChildKey = rightKey;
			largeChildIndex = root * 2 + 1;
		}
		
		// test if root > larger subtree
		if(heap[root] < largeChildKey) {
			// parent < children
			hold = heap[root];
			heap[root] = heap[largeChildIndex];
			heap[largeChildIndex] = hold;
			reheapDown(heap,largeChildIndex, last);
		}
	}
	return;
}

bool insertHeap(int& heap[], int& last, int data) {
	if(last == HEAP_SIZE - 1 )
		return false
	++last;
heap[last] = data;
reheapUp(heap,last);
return true;
}

bool deleteHeap(int heap[], int& last, int& dataOut) {
	if(last < 0)
		return false;
	dataOut = heap[0];
	heap[0] = heap[last];
	last--;
	reheapDown(heap,0,last);
	return true;
}