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

class Graph{
private:
    unsigned int size;
    int* array;
    unsigned int head = 0;
    unsigned int tail = 0;
    unsigned int max;
public:
    Queue(int s);
    ~Queue();
    int front();
    int back();
    void push(int x);
    int pop();
    unsigned int size();
    bool empty();
    bool full();
};

#endif