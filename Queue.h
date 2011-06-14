//  Author: Christopher Lee Jackson


#ifndef __QUEUE_H_INCLUDED__  
#define __QUEUE_H_INCLUDED__  


#include <limits>
#include <iostream>
#include <cmath>


using namespace std;

/*
 *
 *	QUEUE CLASS
 *
 */

class Queue{
private:
    unsigned int size;
    unsigned int head;
    unsigned int tail;
public:
    unsigned int max;
    int* array;
    Queue(int s){array=new int[s]; max=s; size=0; head=0; tail=0;}
    ~Queue(){delete[] array;}
    int front(){return array[head];}
    int back(){return array[tail];}
    void push(int x);
    int pop();
    unsigned int getSize(){return size;}
    bool empty(){return (size==0);}
    bool full(){return (size==max);}
    void reset(){size=0; head=0; tail=0;}
};

inline void Queue::push(int x) {
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
 
inline int Queue::pop() {
    if(empty()) {
        cerr << "Queue is empty. Nothing to pop.";
        exit(-1);
    }
    int temp = head;
    head = (++head)%max;
    size--;
    return array[temp];
}


#endif
