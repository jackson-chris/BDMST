#include <iostream>
#include "BinaryHeap2.h"

using namespace std;

int main() {
    vector<int>* v = new vector<int>();
        v->push_back(2);
        v->push_back(172);
        v->push_back(7432);
        v->push_back(102);
        v->push_back(89);
        v->push_back(383);
        v->push_back(1);
        v->push_back(38394);
        v->push_back(25);
    BinaryHeap* heap = new BinaryHeap(*v);
    while(!heap->isEmpty()) {
        cout << heap->deleteMax() << endl;
    }
    return 0;
}
