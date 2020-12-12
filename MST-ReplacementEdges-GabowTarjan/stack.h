#ifndef _STACK
#define _STACK

#define MAXSTACKSIZE 999999

inline int bitExtract(int num, int a, int b) {
  return (((1 << a)-1) & (num >> (b-1)));
}

int stack[MAXSTACKSIZE];     

int top = -1;            

inline int isempty() {
  return(top == -1);
}
   
inline int isfull() {
  return(top == MAXSTACKSIZE);
}

inline int peek() {
  return stack[top];
}

inline int pop() {
  int data;
	
  if (!isempty()) {
    data = stack[top];
    top = top - 1;   
  } else {
    fprintf(stderr,"Could not retrieve data, Stack is empty.\n");
  }
  return(data);
}

inline void push(int data) {
  if (!isfull()) {
    top = top + 1;   
    stack[top] = data;
#if 1
    fprintf(stdout,"Adding to stack: %9d\n",top);
#endif
  } else {
    fprintf(stderr,"Could not insert data, Stack is full.\n");
  }
  return;
}


#if 0
int bitExtract(int, int, int);
int isempty();
int isfull();
int peek();
int pop();
void push(int);
#endif

#endif
