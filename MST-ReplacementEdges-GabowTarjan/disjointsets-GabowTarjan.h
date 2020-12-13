#ifndef _DISJOINT_GABOW_TARJAN
#define _DISJOINT_GABOW_TARJAN

/* Set the size of a microset */
/* B-1(ceil(log B)) has to be less than a or the word length */
/* and B >= 2 */
#define B  5

typedef struct {
  int Parent;
  int Rank;
} Subset_t;

typedef struct {
  int root;
  int vertices[B];
} MicroSetType;

typedef struct {
  int parent;
  int root;
} MacroSetType;

typedef struct {
  int table[256][8];
} Answer;

typedef struct {
  int *nodetable;
} Table;

typedef struct {
  Table *mark;
} MultiDimTable;

int Find(Subset_t *, int);
void Union(Subset_t *, int, int);
void MicroLink(int*, int*, int, int *, int *);
int find_macromicro(int, int *, int*, MicroSetType *,char answer[4096][32][5], MacroSetType *, int *, int*, int, int*, int*, int*, int);
void makeMacroSet(MicroSetType *, MacroSetType *, int, int*, int);


#endif
