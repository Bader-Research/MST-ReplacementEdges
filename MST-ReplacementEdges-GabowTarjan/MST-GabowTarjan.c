/* Linear Time Algorithm for finding Minimum Spanning Tree Replacement Edges */
/* David A. Bader */
/* July 2019 */

/* Updated December 2020, with Pranhav Sundararajan */
/* Uses Gabow-Tarjan disjoint sets */

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<ctype.h>

#include "disjointsets-GabowTarjan.h"

#define DEBUG_EDGELIST 0
#define DEBUG_ADDEDGE 0
#define DEBUG_GRAPHAL 0
#define DEBUG_BRIDGES 0
#define DEBUG_ROOT 0
#define DEBUG_PATHLABEL 0
#define DEBUG_REPLACE 0
#define DEBUG_MICROSET 0
#define TIMING 0

#define amtMarkTables
char *INFILENAME;
FILE *outfile;
int QUIET;

/***********************************************************************/
#if TIMING
struct timespec timer_start(){
  struct timespec start_time;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start_time);
  return start_time;
}

// call this function to end a timer, returning nanoseconds elapsed as a long
long timer_end(struct timespec start_time){
  struct timespec end_time;
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end_time);
  long diffInNanos = (end_time.tv_sec - start_time.tv_sec) * (long)1e9 + (end_time.tv_nsec - start_time.tv_nsec);
  return diffInNanos;
}
#endif
/***********************************************************************/

void checkPtr(void *ptr) {
  if (ptr == NULL) exit(1);
  return;
}

#if 1
unsigned int nextPowerOfTwo(unsigned int n) {
  unsigned int p;
  p = 1; 
  if (n && !(n & (n - 1))) 
    return n; 
  
  while (p < n)  
    p <<= 1; 
      
  return p; 
}

unsigned int calcPower(unsigned int base, unsigned int exp) {
  int i;
  unsigned int result;

  i = 1;
  result = 1;

  while(i <= exp){
    result *= base;
    i++;
  }

  return result;
}
#endif

/***********************************************************************/
/* Stack */

struct Stack {
  int top;
  int capacity;
  int* array;
};
 
struct Stack* createStack(int capacity) {
  struct Stack* stack = (struct Stack*)malloc(sizeof(struct Stack));
  checkPtr(stack);
  stack->capacity = capacity;
  stack->top = -1;
  stack->array = (int*)malloc(stack->capacity * sizeof(int));
  checkPtr(stack->array);
  return stack;
}
 
void destroyStack(struct Stack* stack) {
  free(stack->array);
  free(stack);
  return;
}
 
int isFull(struct Stack* stack) {
  return stack->top == stack->capacity - 1;
}
 
int isEmpty(struct Stack* stack) {
  return stack->top == -1;
}
 
void push(struct Stack* stack, int item) {
  if (isFull(stack)) {
    fprintf(stderr,"ERROR: Stack is full.\n");
    return;
  }
  stack->array[++stack->top] = item;
  return;
}
 
int pop(struct Stack* stack) {
  if (isEmpty(stack)) {
    fprintf(stderr,"ERROR: Stack is empty.\n");
    return -1;
  }
  return stack->array[stack->top--];
}
 
int peek(struct Stack* stack) {
  if (isEmpty(stack))
    return -1;
  return stack->array[stack->top];
}
 
/***********************************************************************/


/***********************************************************************/
/* MST code */
/***********************************************************************/

typedef struct {
  int v1;
  int v2;
  int w;
  int T;
} Edge_t;

typedef struct {
  int n;
  int m;
  Edge_t* Edgelist;
} Graph_t;

void freeGraph(Graph_t *graph) {
  free(graph->Edgelist);
  free(graph);
  return;
}

Graph_t *CreateGraph(int n, int m) {
  Graph_t *graph;

  graph = (Graph_t *)malloc(sizeof(Graph_t));
  checkPtr(graph);
  graph->n = n;
  graph->m = m;
  graph->Edgelist = (Edge_t *)malloc(m * sizeof(Edge_t));
  checkPtr(graph->Edgelist);

  return graph;
}

int CompareEdges(const void* a, const void* b) {
  return ((Edge_t *)a)->w - ((Edge_t *)b)->w;
}

void PrintMST(Edge_t *result, int e) {
  int i;
  for (i = 0; i < e; i++)
    printf("MST: %6d -- %6d == %6d\n", result[i].v1, result[i].v2, result[i].w);
  return;
}

int getRoot(int *arr, int i) {
  while(arr[i]!=i){
    arr[i]=arr[arr[i]];
    i=arr[i];
  }
  return i;
}

int ConnectedComponents(Graph_t *graph) {
  int c;
  int n, m;
  int x, y, xRoot, yRoot;
  int i;
  int *root;

  n = graph->n;
  m = graph->m;
  c = n;

  root = (int *)malloc (n*sizeof(int));
  checkPtr(root);
  
  for (i=0; i<n; i++)
    root[i]=i;        
  
  for(i=0; i<m; i++){
    x = graph->Edgelist[i].v1;
    y = graph->Edgelist[i].v2;
 
    xRoot = getRoot(root, x);
    yRoot = getRoot(root, y);
 
    if (xRoot!=yRoot) {
      c--;
      root[xRoot]=yRoot;
    }
  }
     
  free(root);
   
  return c;
}

int countLeaves(Graph_t *graph) {
  int n, m;
  int i;
  int leaves;
  int *degrees;

  n = graph->n;
  m = graph->m;

  degrees = (int *)malloc(n*sizeof(int));
  checkPtr(degrees);

  for (i=0; i<n; i++)
    degrees[i]=0;        
  
  for(i=0; i<m; i++){
    degrees[graph->Edgelist[i].v1]++;
    degrees[graph->Edgelist[i].v2]++;
  }

  leaves = 0;
  for (i=0; i<n; i++)
    if (degrees[i]==1) leaves++;

  free(degrees);
   
  return(leaves);
}


void Kruskal(Graph_t *graph) {
  int n;
  Edge_t *result;
  int i;
  int e;
  int x, y, v;
  Subset_t *subsets;
  Edge_t nextEdge;
  
  n = graph->n;
  result = (Edge_t *)malloc(n * sizeof(Edge_t));
  checkPtr(result);

  qsort(graph->Edgelist, graph->m, sizeof(Edge_t), CompareEdges);

  subsets = (Subset_t *) malloc(n * sizeof(Subset_t));
  checkPtr(subsets);

  for (v = 0; v < n; v++) {
      subsets[v].Parent = v;
      subsets[v].Rank = 0;
  }

  i=0; e=0;
  while (e < n - 1) {
    nextEdge = graph->Edgelist[i]; 
    x = Find(subsets, nextEdge.v1);
    y = Find(subsets, nextEdge.v2);
      
    if (x != y) {
      result[e++] = nextEdge;
      Union(subsets, x, y);
      graph->Edgelist[i].T = 1;
    }
    i++;
  }

  free(subsets);
  free(result);
  return;
}




/***********************************************************************/
/* Adjacency List */
/***********************************************************************/

struct node
{
  int vertex;
  int weight;
  struct node* next;
};

typedef struct
{
  int n;
  struct node** adjLists;
} GraphAL_t;


struct node* createNode(int v, int w)
{
  struct node* newNode = (struct node *)malloc(sizeof(struct node));
  checkPtr(newNode);
  newNode->vertex = v;
  newNode->weight = w;
  newNode->next = NULL;
  return newNode;
}
 
GraphAL_t *createGraph_AL(int n)
{
  int i;
  GraphAL_t *graph;
    
  graph = (GraphAL_t *)malloc(sizeof(GraphAL_t));
  checkPtr(graph);
    
  graph->n = n;
 
  graph->adjLists = (struct node**)malloc(n * sizeof(struct node*));
  checkPtr(graph->adjLists);
 
  for (i = 0; i < n; i++)
    graph->adjLists[i] = NULL;
 
  return graph;
}
 
void addEdge(GraphAL_t *graph, int src, int dest, int weight)
{
  struct node* newNode;
  
  /* Add edge from src to dest of weight w */
  newNode = createNode(dest, weight);
  newNode->next = graph->adjLists[src];
  graph->adjLists[src] = newNode;
#if DEBUG_ADDEDGE
  printf("addEdge1: %6d, %6d\n",src, dest);
#endif
  
  /* Add edge from dest to src of weight w */
  newNode = createNode(src, weight);
  newNode->next = graph->adjLists[dest];
  graph->adjLists[dest] = newNode;
  
#if DEBUG_ADDEDGE
  printf("addEdge2: %6d, %6d\n",dest, src);
#endif
  return;
}

void freeGraphAL(GraphAL_t *graph) {
  int v;
  struct node *temp, *nextnode;
      
  for (v = 0; v < graph->n; v++) {
    temp = graph->adjLists[v];
    while (temp) {
      nextnode = temp->next;
      free(temp);
      temp = nextnode;
    }
  }

  free(graph->adjLists);
  free(graph);
  return;
}

 
void printGraphAL(GraphAL_t *graph) {
  int v;
  struct node* temp;
      
  for (v = 0; v < graph->n; v++) {
    temp = graph->adjLists[v];
    printf("printGraphAL: Adjacency list of vertex %6d\n ", v);
    while (temp) {
      printf("%6d -> ", temp->vertex);
      temp = temp->next;
    }
    printf("\n");
  }
  return;
}

/******************************************************/

int min(int a, int b) {
  return ((a<b)?a:b);
}

void dfs_bridges(GraphAL_t *G, int u, int v, int *low, int *pre, int *cnt, int *bridges) {
  int w;
  struct node *p;
  
#if DEBUG_BRIDGES
  printf("dfs_bridges: u:%6d v:%6d\n",u, v);
#endif
  
  pre[v] = *cnt;
  *cnt += 1;
  low[v] = pre[v];
  p=G->adjLists[v];

  while(p!=NULL) {
    w = p->vertex;
    if (pre[w] == -1) {
      dfs_bridges(G, v, w, low, pre, cnt, bridges);
      low[v] = min(low[v], low[w]);
      if (low[w] == pre[w]) {
#if DEBUG_BRIDGES
	printf("%6d - %6d is a bridge\n",v,w);
#endif
	*bridges += 1;
      }
    }
  /* update low number - ignore reverse of edge leading to v */
    else {
      if (w != u)
	low[v] = min(low[v], pre[w]);
    }
    p = p->next;
  }

  return;
}

int countBridges(Graph_t *graph) {
  GraphAL_t *graphAL;
  Edge_t *edges;
  int i, n, m;
  int bridges;

  n = graph->n;
  m = graph->m;
  graphAL = createGraph_AL(n);

  edges = graph->Edgelist;
  for (i=0 ; i<m ; i++)
    addEdge(graphAL, edges[i].v1, edges[i].v2, edges[i].w);

  int cnt;
  int *pre;        /* pre[v] = order in which dfs examines v */
  int *low;        /* low[v] = lowest preorder of any vertex connected to v */

  low = (int *)malloc(n*sizeof(int));
  checkPtr(low);
  pre = (int *)malloc(n*sizeof(int));
  checkPtr(pre);

  for (i = 0; i < n; i++) {
    low[i] = -1;
    pre[i] = -1;
  }

  bridges = 0;

  for (i = 0; i < n; i++)
    if (pre[i] == -1)
      dfs_bridges(graphAL, n-1, n-1, low, pre, &cnt, &bridges);

  free(pre);
  free(low);
  freeGraphAL(graphAL);
  return(bridges);
}

/******************************************************/

/******************************************************/
/* DFS */
/******************************************************/
static int eulerCount;


/******************************************************/
/* Replacement edge PathLabel */
/******************************************************/

#define NULLVERTEX -1

typedef struct {
  int v1;
  int v2;
} Replacement_t;

#define PL_LEFT  0
#define PL_RIGHT 1
#define PL_ANC   2

int PathLabel(int s, int t, int *PARENT,  int *IN, int *OUT,  Replacement_t *REPLACEMENT, Subset_t *sub, int micro[], int number[], MicroSetType *microsub,  char answer[][32][5], MacroSetType *macrosub, int *markTableCounter,  int *amtNodes, int *MicrosetOfroot, int microSetNum, int parent[], int root) {
  int found;
  int k1, k2;
  int v;
  int ehat;
  int PLAN;

#if DEBUG_PATHLABEL
  printf("PathLabel: s:%6d t:%6d\n",s, t);
#endif

  found = 0;

  if ((IN[s] < IN[t]) && (IN[t] < OUT[s])) {
#if DEBUG_PATHLABEL
    printf("PathLabel: PLAN exit\n");
#endif
    return(0);
  }

  if ((IN[t] < IN[s]) && (IN[s] < OUT[t])) {
    PLAN = PL_ANC;
#if DEBUG_PATHLABEL
    printf("PathLabel: PLAN ANC\n");
#endif
    k1 = IN[t];
    k2 = IN[s];
  }
  else {
    if (IN[s] < IN[t]) {
      PLAN = PL_LEFT;
#if DEBUG_PATHLABEL
      printf("PathLabel: PLAN LEFT\n");
#endif
      k1 = OUT[s];
      k2 = IN[t];
    } else {
      PLAN = PL_RIGHT;
#if DEBUG_PATHLABEL
      printf("PathLabel: PLAN RIGHT\n");
#endif
      k1 = OUT[t];
      k2 = IN[s];
    }
  }

#if DEBUG_PATHLABEL
  printf("PathLabel: k1:%6d k2:%6d\n",k1, k2);
#endif

  v = s;
  
#if DEBUG_PATHLABEL
  printf("PathLabel: Tedge: vhat, P[vhat]: <%6d %6d>\n",v, PARENT[v]);
#endif

while (k1 < k2) {
    if (find_macromicro( v, micro, number, microsub,  answer, macrosub, markTableCounter, MicrosetOfroot, microSetNum, parent, amtNodes, PARENT, root) == v) {
    	
		ehat = v;
    	if (REPLACEMENT[ehat].v1 == NULLVERTEX) {			
			REPLACEMENT[ehat].v1 = s;
			REPLACEMENT[ehat].v2 = t;

			found ++;
		}
		

#if DEBUG_REPLACE
	fprintf(outfile,"Replacement edge for <%6d, %6d>: <%6d, %6d>\n",ehat, PARENT[ehat], REPLACEMENT[ehat].v1, REPLACEMENT[ehat].v2);
#endif      

	MicroLink( micro, number,  v, markTableCounter,  amtNodes);
    }

    v = find_macromicro(v, micro, number, microsub, answer, macrosub, markTableCounter, MicrosetOfroot, microSetNum, parent, amtNodes, PARENT, root);

    switch(PLAN) {
    case PL_ANC:   k2 =  IN[v];	break;
    case PL_LEFT:  k1 = OUT[v]; break;
    case PL_RIGHT: k2 =  IN[v]; break;
    }
    
  }


 
  
  return(found);
}

#define MAXLEN 256
#define WRANGE 1024

int checkLine(char *str){
  int c;
  int i;

  c=0;
  i=0;
  while (str[i] != '\0') {
    if (!isspace(str[i]))
      c++;
    i++;
  }

  return c;
}

Graph_t *ReadGraphFromFile(FILE *infile, int *numVertices, int *numEdges) {
  int n, m;
  int s, t, w, zeroBased;
  char sstr[MAXLEN], tstr[MAXLEN], wstr[MAXLEN];
  char lineBuf[MAXLEN];
  int first;
  int i;

  Graph_t *graph;

  zeroBased = 0;
  n=0; m=0;
  first = 1;

  while (!feof(infile)){
    if (fgets(lineBuf, MAXLEN, infile) != NULL) {
      sscanf(lineBuf,"%s %s %s",(char *)&sstr, (char *)&tstr, (char *)&wstr);
      if (checkLine(lineBuf)) {
	if (wstr[0] == 0) {
	  if (!QUIET) {
	    if (first) printf("Using random weights\n");
	  }
	  first = 0;
	}
	s = atoi(sstr);
	t = atoi(tstr);
	if (s>n) n=s;
	if (t>n) n=t;
	if ((s==0)||(t==0)) zeroBased = 1;
	m++;
      }
    }
  }
  n++;

  if (!zeroBased) n--;
  if (!QUIET)
    printf("n: %d m: %d\n",n,m);

  rewind(infile);

  graph = CreateGraph(n, m);

  i = 0;
  while (i<m) {
    fgets(lineBuf, MAXLEN, infile);
    sscanf(lineBuf,"%s %s %s",(char *)&sstr, (char *)&tstr, (char *)&wstr);
    if (checkLine(lineBuf)) {
      s = atoi(sstr);
      t = atoi(tstr);
      if (!first) {
	w = rand() % WRANGE;
      }
      else
	w = atoi(wstr);

      if (!zeroBased) {s--; t--;}
      graph->Edgelist[i].v1 = s;
      graph->Edgelist[i].v2 = t;
      graph->Edgelist[i].w = w;
      graph->Edgelist[i].T = 0;

      i++;
    }
  }
  
  fclose(infile);

  if ((n == 0) || (m == 0)) {
    fprintf(stderr,"ERROR: no edges in graph\n");
    exit(1);
  }
  
  *numVertices = n;
  *numEdges = m;
  return(graph);
}

/******************************************************/

void usage(void) {

  printf("Minimum Spanning Tree Replacement Edges\n\n");
  printf("Usage:\n");
  printf(" -f <filename>   [Input Graph]\n\n");
  printf("Optional arguments:\n");
  /*	printf(" -r              [Use Random Edge Weights]\n"); */
  printf(" -o <filename>   [Output File]\n");
  printf(" -q              [Turn on Quiet mode]\n");
  exit (8);
}

void parseFlags(int argc, char **argv, FILE **infile) {

  if (argc < 2) usage();
  *infile = NULL;
  QUIET = 0;

  while ((argc > 1) && (argv[1][0] == '-')) {

    switch (argv[1][1]) {

    case 'f':
      if (!QUIET)
	printf("Input Graph: %s\n",argv[2]);
      *infile = fopen(argv[2], "r");
      if (*infile == NULL) usage();
      INFILENAME = argv[2];
      argv+=2;
      argc-=2;
      break;

    case 'o':
      if (!QUIET)
	printf("Output file: %s\n",argv[2]);
      outfile = fopen(argv[2], "a");
      if (outfile == NULL) usage();
      argv+=2;
      argc-=2;
      break;

    case 'r':
      if (!QUIET)
	printf("Using random weights\n");
      argv++;
      argc--;
      break;

    case 'q':
      QUIET = 1;
      argv++;
      argc--;
      break;
	
    default:
      fprintf(stderr,"Wrong Argument: %s\n", argv[1]);
      usage();
    }

  }

  if (*infile == NULL) usage();
  
  return;
}

void dfs_microset(GraphAL_t *graph, int v, int * visited,int isRoot,int *d, int b, MicroSetType *microsub, int * microSetNum, int *counter, int micro[], int number[],  int *PARENT,  int *amtNodes, int *MicrosetOfroot, Table *NODE, int *parent, int *IN, int *OUT, struct Stack* stack)
{
    struct node *p;
    int child;
    int i;
    int popnum;
    eulerCount++;
    IN[v] = eulerCount;
    p=graph->adjLists[v];

    visited[v]=1;
#if DEBUG_MICROSET
    printf("Visit current=%d, visited[%d]=%d,isRoot=%d, d[root]=%d,b=%d, microsetnum=%d, counter=%d\n",v, v,visited[v],isRoot,d[v],b, * microSetNum, *counter);
#endif
    while(p!=NULL)
    {
       child=p->vertex;
#if DEBUG_MICROSET
      printf("The current node %d has child %d\n",v,child);
#endif
       if(!visited[child]){
         PARENT[child] = v;
         //fprintf(outfile, "PARENT[%6d]: %6d\n", child,PARENT[child]);
	 dfs_microset(graph, child, visited,-1,d,b,microsub, microSetNum, counter, micro, number,   PARENT,  amtNodes, MicrosetOfroot, NODE, parent, IN, OUT, stack);
       	   if (d[v]<((b+1)/2)) {
#if DEBUG_MICROSET
     printf("Before Updata  d[%d]=%d + d[%d]=%d=%d\n",v,d[v],child,d[child],d[v]+d[child]);
#endif
             d[v]=d[v]+d[child];
	     push(stack, child);

             *counter=*counter+1;
#if DEBUG_MICROSET
     printf("1- Push node %d, into stack as the %d th element, d[%d]=%d, d[%d]=%d\n",child,*counter-1,v,d[v],child,d[child]);
     printf("1 Add node %d, into microset %d, as the %d th element, d[%d]=%d, d[%d]=%d\n",child,*microSetNum,*counter-1,v,d[v],child,d[child]);
#endif
           } else
           {
#if DEBUG_MICROSET
     printf("2- Complete a microset %d, with root  %d and the subtree  has %d elements, d[%d]=%d + d[%d]=%d =%d\n",*microSetNum,v,d[v]+d[child]-2,v,d[v],child,d[child],d[v]+d[child]);
#endif
             popnum=d[v]+d[child]-2;

             NODE[*microSetNum].nodetable = (int *)malloc((popnum+1) * sizeof(int));

             checkPtr(NODE[*microSetNum].nodetable);
            // MARK[*microSetNum].nodetable = (int *)malloc((popnum+1) * sizeof(int));
            MicrosetOfroot[v] = *microSetNum;
            // checkPtr(MARK[*microSetNum].nodetable);
            

             amtNodes[*microSetNum] = popnum;
 			
              i=1;
             for (i=1;i<popnum+1;i++){

	       microsub[*microSetNum].vertices[i]=pop(stack);
                 NODE[*microSetNum].nodetable[i] = microsub[*microSetNum].vertices[i];

                micro[microsub[*microSetNum].vertices[i]] = *microSetNum;
                number[microsub[*microSetNum].vertices[i]] = i;
               // fprintf(outfile, "node[%6d, %6d]: %6d\n",*microSetNum, i,microsub[*microSetNum].vertices[i] );
                //amtNodes[*microSetNum]+=1;
#if DEBUG_MICROSET
     printf("2- The root of microset is %d, the %d th element in microset %d is %d\n",v,i,*microSetNum, microsub[*microSetNum].vertices[i]);
#endif
             }
             
            
             microsub[*microSetNum].root=v;
             
             	for( int j =1; j<amtNodes[*microSetNum]+1;j++)
	{
		
   // fprintf(outfile, "PARENT[%6d]: %6d\n",microsub[*microSetNum].vertices[j], PARENT[microsub[*microSetNum].vertices[j]] );
	 if(micro[PARENT[ microsub[*microSetNum].vertices[j]]]==micro[microsub[*microSetNum].vertices[j]] && (PARENT[microsub[*microSetNum].vertices[j]]!=microsub[*microSetNum].vertices[j]))
	 {

	   parent[*microSetNum] += calcPower(nextPowerOfTwo(b),(j-1)) * number[PARENT[ microsub[*microSetNum].vertices[j]]];
	 //	fprintf(outfile, "parent(%6d,%6d): %6d \n",i,j, parent[i][j]  );  
	 	
	 		
	 }
        
    else{
    	parent[*microSetNum]+=0;
    //	fprintf(outfile, "parent(%6d,%6d): %6d \n",i,j, parent[i][j]  ); 
	} 

        			
	}

             


             *microSetNum=*microSetNum+1;
             *counter=*counter-popnum;
             d[v]=2;
             d[child]=1;
	     push(stack, child);

             *counter=*counter+1;
#if DEBUG_MICROSET
     printf("2- Push node %d, into stack as the %d th element, d[%d]=%d, d[%d]=%d\n",child,*counter-1,v,d[v],child,d[child]);
#endif
           }
       }  
       p=p->next;

    }
    eulerCount++;
    OUT[v] = eulerCount;

       	   if (d[v]>=((b+1)/2)) {
#if DEBUG_MICROSET
    printf("4- Complete a microset %d, with root  %d and the subtree  has %d elements, d[%d]=%d\n",*microSetNum,v,d[v]-1,v,d[v]);
#endif
             popnum=d[v]-1;
             //microsub[*microSetNum].vertices[popnum];

            NODE[*microSetNum].nodetable = (int *)malloc((popnum+1) * sizeof(int));

             checkPtr(NODE[*microSetNum].nodetable);
             
            
             
             
             amtNodes[*microSetNum] = popnum;

          
            i=1;
             for (i=1;i<popnum+1;i++){

	       microsub[*microSetNum].vertices[i]=pop(stack);
                
               
                NODE[*microSetNum].nodetable[i] =  microsub[*microSetNum].vertices[i];
             //   MARK[*microSetNum].nodetable[i] =0;
                micro[microsub[*microSetNum].vertices[i]] = *microSetNum;
                number[microsub[*microSetNum].vertices[i]] = i;
            //amtNodes[*microSetNum]+=1;
             //   fprintf(outfile, "node[%6d, %6d]: %6d\n",*microSetNum, i,microsub[*microSetNum].vertices[i] );
#if DEBUG_MICROSET
     printf("4- The root of microset is %d, the %d th element in microset %d is %d\n",v,i,*microSetNum, microsub[*microSetNum].vertices[i]);
#endif
				
             }
            
            
             microsub[*microSetNum].root=v;
            MicrosetOfroot[v] = *microSetNum;
             for( int j =1; j<amtNodes[*microSetNum]+1;j++)
	{
		
    //fprintf(outfile, "PARENT[%6d]: %6d\n",microsub[*microSetNum].vertices[j], PARENT[microsub[*microSetNum].vertices[j]] );
	 if(micro[PARENT[ microsub[*microSetNum].vertices[j]]]==micro[microsub[*microSetNum].vertices[j]] && (PARENT[ microsub[*microSetNum].vertices[j]]!=microsub[*microSetNum].vertices[j]))
	 {
 
	   parent[*microSetNum] += calcPower(nextPowerOfTwo(b),(j-1)) * number[PARENT[ microsub[*microSetNum].vertices[j]]];
	 //	fprintf(outfile, "parent(%6d,%6d): %6d \n",i,j, parent[i][j]  );  
	 	
	 		
	 }
        
    else{
    	parent[*microSetNum]+=0;
    //	fprintf(outfile, "parent(%6d,%6d): %6d \n",i,j, parent[i][j]  ); 
	} 

        			
	}
             
             *microSetNum=*microSetNum+1;
             *counter=*counter-popnum;
             d[v]=1;

           }


/*
    if (d[v]>=b) {
#if DEBUG_MICROSET
     printf("4 Complete a microset %d, with root  %d and totally has %d elements, d[%d]=%d + d[%d]=%d =%d\n",*microSetNum,v,*counter-1,v,d[v],child,d[child],d[v]+d[child]);
#endif
             popnum=d[v]-1;
             for (i=0;i<popnum;i++){
                microsub[*microSetNum].vertices[i]=pop(&top);
             }
             microsub[*microSetNum].root=v;
             *microSetNum=*microSetNum+1;
             *counter=*counter-popnum;
             d[v]=1;
     }
*/    
    if (isRoot>0) {
             popnum=*counter-1;
        
             NODE[*microSetNum].nodetable = (int *)malloc((popnum+2) * sizeof(int));
             checkPtr(NODE[*microSetNum].nodetable);
          //  MARK[*microSetNum].nodetable = (int *)malloc((popnum+2) * sizeof(int));
           //  checkPtr(MARK[*microSetNum].nodetable);

             //
             amtNodes[*microSetNum] = popnum+1;
            
             for (i=1;i<popnum+1;i++){
             	
                	
	       microsub[*microSetNum].vertices[i]=pop(stack);
                NODE[*microSetNum].nodetable[i] = microsub[*microSetNum].vertices[i];

                micro[microsub[*microSetNum].vertices[i]] = *microSetNum;
                number[microsub[*microSetNum].vertices[i]] = i;
                //amtNodes[*microSetNum]+=1;
               // fprintf(outfile, "node[%6d, %6d]: %6d\n",*microSetNum, i,microsub[*microSetNum].vertices[i] );
#if DEBUG_MICROSET
     printf("5- The root of microset is %d, the %d th element in microset %d is %d\n",-1,i,*microSetNum, microsub[*microSetNum].vertices[i]);
#endif
             }
/* We need to add the root node that is not pushed into the stack*/  
//amtNodes[*microSetNum]+=1;
				microsub[*microSetNum].root=v;  
        MicrosetOfroot[v] = *microSetNum;        
                microsub[*microSetNum].vertices[i]=v;
                NODE[*microSetNum].nodetable[i] = microsub[*microSetNum].vertices[i];
                micro[microsub[*microSetNum].vertices[i]] = *microSetNum;
                number[microsub[*microSetNum].vertices[i]] = i;
               // fprintf(outfile, "node[%6d, %6d]: %6d\n",*microSetNum, i,microsub[*microSetNum].vertices[i] );
                
    for( int j =1; j<amtNodes[*microSetNum]+1;j++)
	{
		
   // fprintf(outfile, "PARENT[%6d]: %6d\n",microsub[*microSetNum].vertices[j], PARENT[microsub[*microSetNum].vertices[j]] );
	 if(micro[PARENT[ microsub[*microSetNum].vertices[j]]]==micro[microsub[*microSetNum].vertices[j]] && (PARENT[ microsub[*microSetNum].vertices[j]]!=microsub[*microSetNum].vertices[j]))
	 {
 
	   parent[*microSetNum] += calcPower(nextPowerOfTwo(b),(j-1)) * number[PARENT[ microsub[*microSetNum].vertices[j]]];
	 //	fprintf(outfile, "parent(%6d,%6d): %6d \n",i,j, parent[i][j]  );  
	 	
	 		
	 }
        
    else{
    	parent[*microSetNum]+=0;
    //	fprintf(outfile, "parent(%6d,%6d): %6d \n",i,j, parent[i][j]  ); 
	} 

        			
	}
				
            
//             microsub[*microSetNum].vertices[popnum]=v;
				
             
#if DEBUG_MICROSET
     printf("5- The root of microset is %d, the %d th element in microset %d is %d\n",-1,i,*microSetNum, microsub[*microSetNum].vertices[i]);
     
     printf("5- Complete a microset %d, with root  %d and totally has %d elements, d[%d]=%d \n",*microSetNum,-1,popnum+1,v,d[v]);
#endif
    }
    //fprintf(outfile, "number[%6d]: %6d\n",8,number[8] );
}
/******************************************************/

inline int bitExtract(int num, int a, int b) {
  return (((1 << a)-1) & (num >> (b-1)));
}


int main(int argc, char **argv) {
  Subset_t *sub;
  int n, m;
  FILE *infile;
  int root;
  int *visited, *PARENT, *IN, *OUT, *L;
  int *NEXT;
  Replacement_t *REPLACEMENT;
  register int i, j, k;
  /*int leaves;*/
  int bridges, maxReplace;
  int cc;
  int replacementFound;
  int *d;
  int b; 
  int  *MicrosetOfroot;
  int microSetNum;
  int numMicro;
  int counter;
  Table *NODE;
  int *amtNodes;
  int *number;
  int *micro;
  int *parent;

  /* 4096 is 2^(ceil(log2 B)*(B-1)) */
  /* 32 is 2^B */
  /* 5 is B */
  char answer[4096][32][5] = {0};

  int c;
  int z;
  struct Stack *stack;

  MicroSetType *microsub;
  MacroSetType *macrosub;

  int *markTableCounter;
  Graph_t *graph;
  Edge_t *edges;
  GraphAL_t *minSpanTree;

  outfile = stdout;
  parseFlags(argc, argv, &infile);

  graph = ReadGraphFromFile(infile, &n, &m);

  cc = ConnectedComponents(graph);
  if (cc>1) {
    fprintf(stderr,"ERROR: Graph has %d Connected Components\n",cc);
    exit(1);
  }

  /*  leaves = countLeaves(graph); */
  bridges = countBridges(graph);

  Kruskal(graph);

#if DEBUG_EDGELIST
  printf("\n");
  for (i=0 ; i<m ; i++) {
    printf("i: %6d, (%6d, %6d) %6d  T:%1d\n",i,graph->Edgelist[i].v1,graph->Edgelist[i].v2, graph->Edgelist[i].w, graph->Edgelist[i].T);
  }
#endif

  minSpanTree = createGraph_AL(n);
  edges = graph->Edgelist;
  for (i=0 ; i<m ; i++) {
    if (graph->Edgelist[i].T) {
#if DEBUG_EDGELIST
      printf("addEdge: <%6d, %6d>\n",edges[i].v1, edges[i].v2);
#endif
      addEdge(minSpanTree, edges[i].v1, edges[i].v2, edges[i].w);
    }
  }
 
  /* Allocate arrays */
  
  visited = (int *)malloc(n * sizeof(int));
  checkPtr(visited);

  PARENT = (int *)malloc(n * sizeof(int));
  checkPtr(PARENT);

  IN = (int *)malloc(n * sizeof(int));
  checkPtr(IN);
  
  sub = (Subset_t *)malloc(n * sizeof(Subset_t));
  checkPtr(sub);
  
  MicrosetOfroot = (int *)malloc(n * sizeof(int));
  checkPtr(MicrosetOfroot);
  
  microsub = (MicroSetType *)malloc((n+1) * sizeof(MicroSetType));
  checkPtr(microsub);

  d = (int *)malloc(n * sizeof(int));
  checkPtr(d);

  OUT = (int *)malloc(n * sizeof(int));
  checkPtr(OUT);

  NEXT = (int *)malloc(n * sizeof(int));
  checkPtr(NEXT);
  
  L = (int *)malloc(n * sizeof(int));
  checkPtr(L);

  REPLACEMENT = (Replacement_t *)malloc(n * sizeof(Replacement_t));
  checkPtr(REPLACEMENT);

  number = (int *)malloc((n+1) * sizeof(int));
  checkPtr(number);

  micro = (int *)malloc((n+1) * sizeof(int));
  checkPtr(micro);

#if TIMING
  int REPEAT = (1<<23) / n;
  if (REPEAT < 10) REPEAT = 10;
  if (!QUIET)
    printf("REP %d\n",REPEAT);
  struct timespec rtime = timer_start();
  for (r=0; r<1 ; r++){
#endif

    /* pick random vertex as root */
    srand(time(0));
    root = rand() % n;
#if DEBUG_ROOT
    printf("\n");
    printf("Root vertex: %6d\n",root);
#endif

    /* root tree */
    /* perform DFS at root, setting Parent(v), IN(v), OUT(v) */

    for (i=0 ; i<n ; i++)
      visited[i] = 0;

    PARENT[root] = root;
    eulerCount = 0;
     
    for (i=0 ; i<n ; i++)
      d[i] = 1;
   
    b = B;
    numMicro = 2*n/(b-1) + 1;  
  
    NODE = (Table *)malloc((numMicro+1) * sizeof(Table));
    checkPtr(NODE);
   
    amtNodes = (int *)malloc((numMicro+1)* sizeof(int));
    checkPtr(amtNodes);
  
    markTableCounter = (int *)malloc((numMicro+1) * sizeof(int));
    checkPtr(markTableCounter);

    parent = (int *)malloc((numMicro+1) * sizeof(int));
    checkPtr(parent);

    microSetNum = 1;
    counter = 1;

    stack = createStack(n);

    for (i=0 ; i<n ; i++) {
      visited[i] = 0;
      micro[i] = -1;
      number[i]=0;
      MicrosetOfroot[i] = -1;
    }

    for (i=1 ; i<numMicro+1 ; i++) {
      parent[i] = 0;
      amtNodes[i] = 0;
      markTableCounter[i] = 0;
    }

    dfs_microset(minSpanTree, root,visited, 1,d, b, microsub, & microSetNum, &counter, micro, number,  PARENT,  amtNodes, MicrosetOfroot, NODE, parent, IN, OUT, stack);

    macrosub = (MacroSetType *)malloc((microSetNum+1) * sizeof(MacroSetType));
    checkPtr(macrosub);

    j=0; 
    PARENT[root] = root;

    makeMacroSet(microsub, macrosub, (microSetNum+1), micro, root);

#if 0
    const int parameter = numMicro+1;
    const int numMarkTables = (int)pow(2,parameter);
    const int numParentTables = (int)pow(2,(int)((parameter-1)*(int)ceil(log2(b))));
#endif

    /* initialize answer table */
    z = 1;

    for (i=1 ; i<microSetNum+1 ; i++) {
      for(j=0; j<32 ; j++) {
	for(k=1 ; k<b ; k++) {
	  c = j >>(k-1) & 1;
	  if(c==0)
	    answer[parent[i]][j][k] = k;
	  else {
	    z=k;
	    c=1;
	    while (z!=0 && c!=0) {
	      z = bitExtract(parent[i],3,(z)*3-2);
	      if (z!=0)
		c = j >>(z-1) & 1;
	      else
		break;
	    }
	    if (z==0)
	      answer[parent[i]][j][k] = 0;
	    if(c==0 && z!=0)
	      answer[parent[i]][j][k] = z;
	  }
	}
      }
    }


  /**********************************************************/
  /* Replacement Edges */
  /**********************************************************/

  /* Create solution REPLACEMENT[], where REPLACEMENT[i] is the replacement edge for MST edge REPLACEMENT[i], REPLACEMENT[PARENT[i]] */
    for (i=0 ; i<n ; i++) {
      REPLACEMENT[i].v1 = NULLVERTEX;
      /* Lazy: no need to init v2 to NULL */
      /* REPLACEMENT[i].v2 = NULLVERTEX; */
    }
  
  /* Scan through m-n+1 non-tree edges, and run PathLabel for each */
    maxReplace = (n-1) - bridges; 
    replacementFound = 0;

    for (i=0 ; (i<m) && (replacementFound < maxReplace) ; i++) {
      if (edges[i].T == 0) {
#if DEBUG_PATHLABEL
	printf("PathLabel for edge <%6d, %6d>\n", edges[i].v1, edges[i].v2);
#endif
	replacementFound += PathLabel(edges[i].v1, edges[i].v2, PARENT,  IN, OUT,  REPLACEMENT, sub,   micro, number, microsub,  answer, macrosub, markTableCounter, amtNodes, MicrosetOfroot, microSetNum, parent, root);
	replacementFound += PathLabel(edges[i].v2, edges[i].v1, PARENT,  IN, OUT,  REPLACEMENT, sub,   micro, number, microsub,   answer, macrosub, markTableCounter,  amtNodes, MicrosetOfroot, microSetNum, parent,root);
      }
    }

 
#if TIMING
  }
  long time_elapsed_nanos = timer_end(rtime);
  fprintf(outfile,"%s %d %d %9.6f\n",INFILENAME, n,m, (double)time_elapsed_nanos / (double)REPEAT);
#endif
#if DEBUG_GRAPHAL
  printGraphAL(minSpanTree);
#endif

 if (!QUIET) {
   for (i=0 ; i<n ; i++) {
     if (REPLACEMENT[i].v1 != NULLVERTEX)
       fprintf(outfile,"Replacement of MST edge <%6d, %6d>: <%6d, %6d>\n",i, PARENT[i], REPLACEMENT[i].v1, REPLACEMENT[i].v2);
   }
 } 
 fclose(outfile);

 destroyStack(stack);
 free(graph);
 freeGraphAL(minSpanTree);
 free(REPLACEMENT); 
 free(OUT);
 free(IN);
 free(PARENT);
 free(visited);
 free(sub);
 free(NODE);
 free(parent);
 free(macrosub);
 free(microsub);
 free(MicrosetOfroot);
 free(markTableCounter);
 free(amtNodes);
 free(micro);
 free(number);
 return 0;
}

