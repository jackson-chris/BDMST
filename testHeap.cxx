#include <iostream>
#include "BinaryHeap2.h"

using namespace std;

int main() {
    vector<int>* v = new vector<int>();
    for(int n = 10; n >= 0; n--)
        v->pushBack(n);
    BinaryHeap* heap = new BinaryHeap(v);
    while(!heap->isEmpty()) {
        cout << heap->deleteMax();
    }
    return 0;
}