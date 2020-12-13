#include<stdlib.h>
#include "stack.h"

int MAXSTACKSIZE;
int *stack;

int stack_init(int n, int *top) {
  MAXSTACKSIZE = 999*n;
  stack = (int *)(MAXSTACKSIZE * sizeof(int));
  *top = -1;
  return (stack != NULL);
}

void stack_free() {
  free(stack);
  return;
}

