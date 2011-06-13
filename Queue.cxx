#include "Queue.h"

Queue::Queue(int s) {
    array = new int[s];
    max = s;
    size = 0;
}

Queue::~Queue() {
    delete array;
}

int Queue::front() {
    return array[head];
}

int Queue::back() {
    return array[tail];
}

void Queue::push(int x) {
    if(empty()) {
        size++;
        array[tail] = x;
    }
    else if (full()){
        // size stays the same but the tail shifts along with the head.
        array[(++tail)%max] = x;
        head = (++head)%max;
    }
    else {
        size++;
        array[(++tail)%max] = x;
    }
}

int Queue::pop() {
    if(empty()) {
        cerr << "Queue is empty. Nothing to pop.";
        exit -1;
    }
    int temp = head;
    head = (++head)%max;
    size--;
    return array[temp];
}

unsigned int Queue::size() {
    return size;
}

bool Queue::empty() {
    return ( size == 0 );
}

bool Queue::full() {
    return ( size == max );
}