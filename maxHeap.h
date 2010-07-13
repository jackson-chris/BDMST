//  Author: Christopher Lee Jackson
//  

#include <limits>
#include <iostream>
#include <vector>
#include "Graph.h"







template <typename Comparable>
class BinaryHeap {
public:
	explicit BinaryHeap (int capacity = 100);
	explicit BinaryHeap (const vector<Comparable> & items);
	
	bool isEmpty() const;
	const Comparable & findMin() const;
	
	void insert( const Comparable & x);
	void deleteMax();
	void deleteMax( Comparable & maxItem);
	void makeEmpty();
private:
	int currentSize; // Number of elements in heap
	vector<Comparable> array; // the heap array
	void buildHeap();
	void percolateDown(int hole);
};
