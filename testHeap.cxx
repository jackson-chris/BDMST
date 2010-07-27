#include <iostream>
#include "BinaryHeap.h"
#include "Graph.h"

using namespace std;

int main() {
    BinaryHeap* heap = new BinaryHeap(false);
    int i = 0;
    while (i < 4) {
        Hub* d = new Hub();
        d->vertId = i++;
        for(int x = 0; x < i; x++) {
            d->edges.push_back(new Edge());
        }
        heap->insert(d);
    }
    while(!heap->isEmpty()) {
        cout << heap->deleteMax()->edges.size() << endl;
    }
    return 0;
}
