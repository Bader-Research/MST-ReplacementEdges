#include "disjointsets-GabowTarjan.h"
#include<stdio.h>
#include<math.h>

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


void MicroLink(int micro[], int number[], int v, int *markTableCounter, int * amtNodes) {
  int z;
#if 0
  int x, y;
#endif
  
#if 0
  x = amtNodes[micro[v]];
  y = x-number[v];
#endif
  z = number[v];
  {
    markTableCounter[micro[v]] += (int) (pow(2,(z-1)) +0.5);
  }
  return;
}

int microfind(int v, int micro[], int number[], MicroSetType *microsub, char answer[4096][32][5], int *markTableCounter, int parent[], int *PARENT ) {
  int i,j,k;
  i = micro[v];
  j = number[v];
  k = answer[parent[i]][markTableCounter[i]][j];

  if (k==0) {
    if(micro[PARENT[v]]==i || PARENT[v]==microsub[i].root)
      return microsub[i].root;
    else
      return PARENT[v];
  }
  else
    return microsub[i].vertices[k];
}

int macrofind(int x, MacroSetType *macrosub, int micro[], int MicrosetOfroot[],MicroSetType *microsub, int PARENT[], int root, int number[], int *markTableCounter, int amtNodes[], char answer[][32][5], int *parent ) {
  if (microfind(x, micro, number, microsub,  answer,markTableCounter,parent,PARENT )!=x) {
    if (micro[PARENT[x]]!=micro[x] || PARENT[x]!=microsub[micro[x]].root)
      return PARENT[x];
  }

  /* DAB COMMENT: Not needed? */
#if 0
  if(macrosub[MicrosetOfroot[x]].parent == macrosub[MicrosetOfroot[x]].root) {
    return macrosub[MicrosetOfroot[x]].parent;
  }
#endif  
  if(macrosub[MicrosetOfroot[x]].parent != macrosub[MicrosetOfroot[x]].root) {
    macrosub[MicrosetOfroot[x]].parent=macrofind(macrosub[MicrosetOfroot[x]].parent,macrosub,micro, MicrosetOfroot,microsub, PARENT,root,number, markTableCounter, amtNodes,answer, parent  );
  }
  return macrosub[MicrosetOfroot[x]].parent;
}

void macrounite(int x, int y, MacroSetType *macrosub, int micro[], int MicrosetOfroot[], MicroSetType *microsub, int *PARENT, int root, int *number, int*markTableCounter, int*amtNodes, char answer[4096][32][5], int *parent) {
  x = PARENT[x];
  while (MicrosetOfroot[x]==-1) {
    x = PARENT[x];
  }
  macrosub[MicrosetOfroot[y]].parent =macrosub[MicrosetOfroot[x]].parent;
  return;
}

void makeMacroSet(MicroSetType *microsub, MacroSetType *macrosub, int n, int micro[], int root) {
  int i;
  for(i=1; i<n; i++) {
    macrosub[i].parent = microsub[i].root;
    macrosub[i].root = microsub[i].root;
  }
}

int find_macromicro(int v, int *micro, int number[], MicroSetType *microsub,  char answer[4096][32][5], MacroSetType *macrosub, int *markTableCounter, int MicrosetOfroot[],int microSetNum, int parent[] , int amtNodes[], int PARENT[], int root) {
  int x, y;
  x=v;
  if(micro[x]!=micro[microfind(x, micro, number, microsub,  answer,markTableCounter,parent,PARENT )]) {
    y = microfind(x, micro, number, microsub,  answer,markTableCounter, parent,PARENT);
    if (MicrosetOfroot[y]==-1)
      return y;
    x= macrofind(y,macrosub, micro, MicrosetOfroot,microsub,PARENT,root,number, markTableCounter, amtNodes,answer, parent);
    if (MicrosetOfroot[x]==-1)
      return x;
    while (micro[x]!=micro[microfind(x, micro, number, microsub, answer,markTableCounter, parent,PARENT)]) {
      macrounite(macrofind(x,macrosub,micro,MicrosetOfroot,microsub, PARENT, root,number, markTableCounter, amtNodes, answer, parent),x, macrosub, micro, MicrosetOfroot,microsub, PARENT,root, number, markTableCounter, amtNodes, answer,parent);
      if(MicrosetOfroot[x]==-1)
        return x;
      x = macrofind(x,macrosub,micro,MicrosetOfroot,microsub, PARENT, root,number, markTableCounter, amtNodes, answer, parent);
      if(MicrosetOfroot[x]==-1)
        return x;
    }
  }
  return microfind(x, micro, number, microsub, answer,markTableCounter, parent,PARENT);
}

void makeSet(int n, Subset_t *sub) {
  int i;
  for(i=0; i<n; i++) {
    sub[i].Parent =  i;
    sub[i].Rank = 0;
  }
  return;
}


