#ifndef _STACK
#define _STACK

#include<stdio.h>

#define INLINE inline

extern int MAXSTACKSIZE;
extern int *stack;


int stack_init(int, int*);
void stack_free();

INLINE int isempty(int top) {
  return(top == -1);
}
   
INLINE int isfull(int top) {
  return(top == MAXSTACKSIZE - 1);
}

INLINE int peek(int top) {
  return stack[top];
}

INLINE int pop(int *top) {
  int data;
	
  if (!isempty(*top)) {
    data = stack[*top];
    *top = *top - 1;   
  } else {
    fprintf(stderr,"Could not retrieve data, Stack is empty.\n");
  }
  return(data);
}

INLINE void push(int data, int *top) {
  if (!isfull(*top)) {
    *top = *top + 1;   
    stack[*top] = data;
  } else {
    fprintf(stderr,"Could not insert data, Stack is full.\n");
  }
  return;
}

#endif
