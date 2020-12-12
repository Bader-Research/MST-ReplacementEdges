#include "disjointsets-Tarjan.h"

int Find(Subset_t *subsets, int i) {
  if (subsets[i].Parent != i)
    subsets[i].Parent = Find(subsets, subsets[i].Parent);

  return subsets[i].Parent;
}

void Union(Subset_t *subsets, int x, int y) {
  register int xroot, yroot;

  Subset_t *sx, *sy;
  
  xroot = Find(subsets, x);
  yroot = Find(subsets, y);

  sx = subsets+xroot;
  sy = subsets+yroot;
  
  if (sx->Rank < sy->Rank)
    sx->Parent = yroot;
  else if (sx->Rank > sy->Rank)
    sy->Parent = xroot;
  else {
    sy->Parent = xroot;
    sx->Rank++;
  }
  return;
}

void Link(Subset_t *subsets, int x, int y) {
  register int xroot, yroot; 
  xroot = Find(subsets, x);
  yroot = Find(subsets, y); 

  subsets[y].Parent = xroot;
  subsets[x].Rank++;
  
  return;
}

void makeSet(int n, Subset_t *sub) {
  int i;
  for(i=0; i<n; i++) {
    sub[i].Parent =  i;
    sub[i].Rank = 0;
  }
  return;
}


