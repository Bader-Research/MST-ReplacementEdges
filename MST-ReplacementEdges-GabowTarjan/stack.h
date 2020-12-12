#ifndef _STACK
#define _STACK

#include<stdio.h>

int MAXSTACKSIZE;
int *stack;
int stack_init(int, int*);
void stack_free();

inline int bitExtract(int num, int a, int b) {
  return (((1 << a)-1) & (num >> (b-1)));
}

inline int isempty(int top) {
  return(top == -1);
}
   
inline int isfull(int top) {
  return(top == MAXSTACKSIZE - 1);
}

inline int peek(int top) {
  return stack[top];
}

inline int pop(int *top) {
  int data;
	
  if (!isempty(*top)) {
    data = stack[*top];
    *top = *top - 1;   
  } else {
    fprintf(stderr,"Could not retrieve data, Stack is empty.\n");
  }
  return(data);
}

inline void push(int data, int *top) {
  if (!isfull(*top)) {
    *top = *top + 1;   
    stack[*top] = data;
  } else {
    fprintf(stderr,"Could not insert data, Stack is full.\n");
  }
  return;
}

#endif
