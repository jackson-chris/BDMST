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
    if(full()) {
        cerr << "Queue is full.";
        exit -1;
    }
    if(empty()) {
        size++;
        array[tail] = x;
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