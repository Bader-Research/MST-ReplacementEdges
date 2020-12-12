#ifndef _DISJOINT_TARJAN
#define _DISJOINT_TARJAN

typedef struct {
  int Parent;
  int Rank;
} Subset_t;

int Find(Subset_t *, int);
void Union(Subset_t *, int, int);
void Link(Subset_t *, int, int);
void makeSet(int, Subset_t *);

#endif
