/**
 * snarkhunter.c
 *
 * Snarkhunter: a generator for cubic graphs and snarks.
 *
 * The latest version of snarkhunter can be found here:
 * http://caagt.ugent.be/cubic/
 *
 * Author: Jan Goedgebeur (jan.goedgebeur@ugent.be)
 * In collaboration with Gunnar Brinkmann and Brendan McKay
 *
 */

#ifndef _SNARKHUNTER_H
#define	_SNARKHUNTER_H

/********************************DEFINES***************************************/

/******************************Configuration***********************************/

//Uncomment to turn debugmode on
//#define _DEBUG


/**
 * Before outputting graphs, they are first written to a buffer. These macros
 * determine the size of those buffers.
 * For snarks a different buffersize is used (which should be much smaller).
 */
#define LISTLENGTH 10000
#define LISTLENGTH_SNARKS 100


/***************You should not change anything below this line*****************/


/*****************************Some useful macros*******************************/

#if WORDSIZE == 64
#define BIT(i) (1ULL << (i))
#else //i.e. WORDSIZE == 32
#define BIT(i) (1U << (i))
#endif

//Max 42 vertices, since we also use edge labels (42*3/2=63)
#define MAXN_EXTRA 42

#define BIT64(i) (1ULL << (i))

#define MAX(a, b) (((a) >= (b)) ? (a) : (b))
#define MIN(a, b) (((a) <= (b)) ? (a) : (b))


//Uncomment to use builtin popcount (else uses nauty popcount)
//Need to compile with gcc -mpopcnt or -msse4.2 to enable __builtin_popcount!
//USE_BUILTIN_POPCOUNT is significantly faster for 64bit
#define USE_BUILTIN_POPCOUNT

#ifdef USE_BUILTIN_POPCOUNT

#if WORDSIZE == 64
#define POPC(x) (__builtin_popcountll(x)) 
#else //i.e. WORDSIZE == 32
#define POPC(x) (__builtin_popcount(x)) 
#endif

#define POPCLL(x) (__builtin_popcountll(x)) 

#else
#define POPC(x) POPCOUNT(x) /* needs nauty.h */
#endif

//Macros for debugging
#ifdef _DEBUG

#define DEBUGMSG(msg) fprintf(stderr, "%s\n", msg);

#define DEBUGASSERT(assertion) if(!(assertion)) {fprintf(stderr, "%s:%u Assertion failed: %s\n", __FILE__, __LINE__, #assertion); exit(1);}

#else

#define DEBUGMSG(msg)

#define DEBUGASSERT(assertion)

#endif

/*******************************Other defines**********************************/

#define REG 3

/**
 * Priority of operations
 * ----------------------
 * Generation of prime graphs:
 * 1. lollipop_diamond
 * 2. nonadj_edge_diamond
 * 3. edge_diamond
 *
 * Generation of (remaining) cubic graphs:
 * 1. triangle
 * 2. nonadj_edge
 */

#define EDGE_INSERTED 1
#define TRIANGLE_INSERTED 0
#define INIT -1
#define MAJOR_EDGE_INSERTED 2

#define EDGE_DIAMOND_INSERTED 3 /* For the generation of the irreducible graphs */
#define NONADJ_EDGE_DIAMOND_INSERTED 4 /* For the generation of the irreducible graphs */
#define LOLLIPOP_DIAMOND_INSERTED 5 /* For the generation of the irreducible graphs */

/* Nauty worksize */
#define WORKSIZE 50 * MAXM


//The max and min of the various (nauty) colour
#define MAX_EDGE_COLOUR_TWO 14

#define MAX_EDGE_COLOUR_DISTANCE_THREE 30
#define MIN_EDGE_COLOUR_DISTANCE_THREE 6

#define MAX_EDGE_COLOUR_COMBINATION (4*MAX_EDGE_COLOUR_TWO + 1 + 4*MAX_EDGE_COLOUR_DISTANCE_THREE + 1)

//Vertex colours for the partitions
#define MAX_VERTEX_COLOUR_DISTANCE_TWO (3 * MAX_EDGE_COLOUR_TWO)
#define MAX_VERTEX_COLOUR_DISTANCE_THREE (3 * MAX_EDGE_COLOUR_DISTANCE_THREE)

//3*7+1
#define MAX_VERTEX_COLOUR_DISTANCE_THREE_TRIPOD 22

//3*15+1
#define MAX_VERTEX_COLOUR_DISTANCE_FOUR_TRIPOD 46

//For the timing
#define time_factor sysconf(_SC_CLK_TCK)


/* Some typedefs */
//Last element is not used. Is about 1% faster than using GRAPH[MAXN][REG].
typedef unsigned char GRAPH[MAXN][REG + 1];

//2 edges:  (0,1) and (2,3)
typedef unsigned char EDGEPAIR[4];
typedef unsigned char SQUARE[4];
typedef unsigned char TRIANGLE[4]; //4 instead of 3, because it's slightly more efficient
typedef unsigned char EDGE[2];

typedef unsigned char EDGETRIPLE[6];

typedef unsigned char EDGE4TUPLE[8];

/* An irred triangle always has the shape of a diamon and consists out of 2 triangles */
typedef unsigned char IRRED_TRIANGLE[4];


/****************************Global variables**********************************/

unsigned long long int number_of_graphs_found = 0;
int number_of_vertices;
int min_order = 0; /* The min order of a graph with the specified girth */

GRAPH current_graph;
int current_number_of_vertices;
int current_number_of_edges;

EDGE edgelist[(3 * MAXN) / 2];
int edgelist_size = 0;

/**
 * The index of each edge in the edgelist.
 * Is used to determine the corresponding last edge more efficiently.
 */
unsigned char edge_index[MAXN][MAXN];

//A mapping from each edge to a unique label (used in determine_edgepair_orbits)
unsigned char edge_labels[MAXN][MAXN];

//Uses edge_labels as index for the rows and columns
//Can't use unsigned char, because there may be more than 255 edgepairs
unsigned int edgepair_index[3*MAXN/2][3*MAXN/2];


unsigned int edgetriple_index[3*MAXN/2][3*MAXN/2][3*MAXN/2];

unsigned int edgetriple_index_g7[3*MAXN/2][3*MAXN/2][3*MAXN/2];

/**
 * To determine the index of a vertexset in the list of vertexsets efficiently.
 * Using malloc because it's safer (can exit if not enough memory can be allocated)
 */
//unsigned int vertexset_index[MAXN][MAXN][MAXN][MAXN];
unsigned int (*vertexset_index)[MAXN][MAXN][MAXN][MAXN];

//Using malloc since this is safer!
unsigned int (*edge4tuple_index)[3*MAXN/2][3*MAXN/2][3*MAXN/2][3*MAXN/2];

//Using malloc since this is safer!
unsigned int (*edge4tuple_index_g7)[3*MAXN/2][3*MAXN/2][3*MAXN/2][3*MAXN/2];

//The edge by which the inserted edge was rejected in the previous iteration
EDGE previous_min_edge;
EDGE previous_min_edge_g7;
EDGE min_edge;
int min_edge_is_part_of_square = 0;

unsigned char previous_min_vertex = 0;

//The edgepairs with index < edgepair_list_index_squares will only yield graphs with girth >= 5
int edgepair_list_index_squares = 0;

EDGE min_edges[(3 * MAXN) / 2];
int min_edges_size = 0;
setword min_edge_bitvectors[(3 * MAXN) / 2];

unsigned char min_vertices_tripod[MAXN];
int min_vertices_tripod_size = 0;

//The list of reducible triangles (i.e. a triangle with three different neighbours)
TRIANGLE reducible_triangles[MAXN];
int number_of_reducible_triangles = 0;

//The list of irreducible triangles (a.k.a. diamonds) of the current graph (i.e. a triangle with at least one common neighbour)
IRRED_TRIANGLE irreducible_triangles[MAXN];
int number_of_irreducible_triangles = 0;
setword irreducible_triangles_bitvector = 0;

EDGE bridges[MAXN]; //Better bound possible
int number_of_bridges = 0;

//is_bridge[i][j] = 1 if (i,j) is a bridge, else value is 0
//Wanring: i is supposed to be < j!
unsigned char is_bridge[MAXN][MAXN];

//Is 1 if all edges of current_graph are reducible
unsigned char all_edges_are_reducible = 0;

//Is only valid for tripod operation!
unsigned char cutvertices[MAXN];
int number_of_cutvertices = 0;

unsigned char is_cutvertex[MAXN];

int number_of_vertices_marked_cutvertex = 0;

int only_output_2_connected_graphs = 0;

int parent_is_3connected = 0;


//Some variables for the find_squares-method
SQUARE squares_global[4];
int squares_global_size = 0;

setword squares_global_bitvectors[4];

EDGE adjacent_squares_global[2];
int adjacent_squares_global_size = 0;


/* Variables needed for the generation of the irreducible graphs */
int max_number_of_vertices_irred_graph = 0;

IRRED_TRIANGLE edge_diamonds[MAXN];
int number_of_edge_diamonds = 0;

IRRED_TRIANGLE nonadj_edge_diamonds[MAXN];
int number_of_nonadj_edge_diamonds = 0;

IRRED_TRIANGLE lollipop_diamonds[MAXN];
int number_of_lollipop_diamonds = 0;

//Eligible edges are edges that aren't fully in a diamond
EDGE eligible_edges[(3 * MAXN) / 2];
int eligible_edges_size = 0;


//Variables only used in case of snarks on level n - 2: here the calls to nauty are delayed
int lab_parent[MAXN];
int ptn_parent[MAXN];

unsigned char colours_three_parent[MAXN][MAXN];
unsigned char min_colour_three_parent;

EDGE edgelist_parent[(3 * MAXN) / 2];
int edgelist_parent_size;

unsigned char edge_index_parent[MAXN][MAXN];


/* Variables for the (nauty) colours */
unsigned char colours_one[MAXN][MAXN];
unsigned char min_colour_one;

unsigned char colours_two[MAXN][MAXN];
unsigned char min_colour_two;

unsigned char colours_three[MAXN][MAXN];
unsigned char min_colour_three;

//Using colours_two for now!
unsigned char colours_h_operation[MAXN][MAXN];
unsigned char min_colour_h_operation;

//unsigned char colours_four[MAXN][MAXN]; //storing colours_four doesn't really help
//unsigned char min_colour_four;

//Used in has_min_colour_distance_three
unsigned char frequencies[MAX_EDGE_COLOUR_DISTANCE_THREE + 1];


setword vertex_colours_long_two[MAXN];
setword vertex_colours_long_three[MAXN];
setword vertex_colours_long_four[MAXN];

unsigned char vertex_colours_long_three_popc[MAXN];
unsigned char vertex_colours_long_four_popc[MAXN];
unsigned char min_colour_four_popc;

setword vertex_neighbourhood[MAXN];

//arrays for determine_vertex_partitions
//Note that +1 is necessary here
unsigned char vertex_colours_one[MAX_VERTEX_COLOUR_DISTANCE_TWO + 1][MAXN];
unsigned char vertex_colours_one_size[MAX_VERTEX_COLOUR_DISTANCE_TWO + 1];

unsigned char vertex_colours_marked[MAX_VERTEX_COLOUR_DISTANCE_THREE + 1][MAXN];
unsigned char vertex_colours_marked_size[MAX_VERTEX_COLOUR_DISTANCE_THREE + 1];

//arrays for determine_vertex_partitions
//Note that +1 is necessary here
unsigned char vertex_colours_three_tripod[MAX_VERTEX_COLOUR_DISTANCE_THREE_TRIPOD + 1][MAXN];
unsigned char vertex_colours_three_tripod_size[MAX_VERTEX_COLOUR_DISTANCE_THREE_TRIPOD + 1];

unsigned char vertex_colours_four_tripod[MAX_VERTEX_COLOUR_DISTANCE_FOUR_TRIPOD + 1][MAXN];
unsigned char vertex_colours_four_tripod_size[MAX_VERTEX_COLOUR_DISTANCE_FOUR_TRIPOD + 1];


/* Some marks */
unsigned int marks[MAXN];

#define MAXVAL INT_MAX - 1
static int markvalue = MAXVAL;
#define RESETMARKS {int mki; if ((markvalue += 1) > MAXVAL) \
      { markvalue = 1; for (mki=0;mki<MAXN;++mki) marks[mki]=0;}}
#define MARK(v) marks[v] = markvalue
#define UNMARK(v) marks[v] = markvalue - 1
#define ISMARKED(v) (marks[v] == markvalue)

unsigned int marks2[MAXN];
static int markvalue2 = MAXVAL;
#define RESETMARKS2 {int mki; if ((markvalue2 += 1) > MAXVAL) \
      { markvalue2 = 1; for (mki=0;mki<MAXN;++mki) marks2[mki]=0;}}
#define MARK2(v) marks2[v] = markvalue2
#define UNMARK2(v) marks2[v] = markvalue2 - 1
#define ISMARKED2(v) (marks2[v] == markvalue2)

unsigned int marks3[MAXN];
static int markvalue3 = MAXVAL;
#define RESETMARKS3 {int mki; if ((markvalue3 += 1) > MAXVAL) \
      { markvalue3 = 1; for (mki=0;mki<MAXN;++mki) marks3[mki]=0;}}
#define MARK3(v) marks3[v] = markvalue3
#define UNMARK3(v) marks3[v] = markvalue3 - 1
#define ISMARKED3(v) (marks3[v] == markvalue3)

//Marks to determine the nauty partitions
unsigned int marks_vertex_colour_2[MAX_VERTEX_COLOUR_DISTANCE_TWO + 1];
static int markvalue_vertex_colour_2 = MAXVAL;
#define RESETMARKS_VERTEX_COLOUR_2 {int mki; if ((markvalue_vertex_colour_2 += 1) > MAXVAL) \
      { markvalue_vertex_colour_2 = 1; for (mki=0;mki<MAX_VERTEX_COLOUR_DISTANCE_TWO + 1;++mki) marks_vertex_colour_2[mki]=0;}}
#define MARK_VERTEX_COLOUR_2(v) marks_vertex_colour_2[v] = markvalue_vertex_colour_2
#define ISMARKED_VERTEX_COLOUR_2(v) (marks_vertex_colour_2[v] == markvalue_vertex_colour_2)


unsigned int marks_vertex_colour_3[MAX_VERTEX_COLOUR_DISTANCE_THREE + 1];
static int markvalue_vertex_colour_3 = MAXVAL;
#define RESETMARKS_VERTEX_COLOUR_3 {int mki; if ((markvalue_vertex_colour_3 += 1) > MAXVAL) \
      { markvalue_vertex_colour_3 = 1; for (mki=0;mki<MAX_VERTEX_COLOUR_DISTANCE_THREE + 1;++mki) marks_vertex_colour_3[mki]=0;}}
#define MARK_VERTEX_COLOUR_3(v) marks_vertex_colour_3[v] = markvalue_vertex_colour_3
#define ISMARKED_VERTEX_COLOUR_3(v) (marks_vertex_colour_3[v] == markvalue_vertex_colour_3)


unsigned int marks_vertex_colour_three_tripod[MAX_VERTEX_COLOUR_DISTANCE_THREE_TRIPOD + 1];
static int markvalue_vertex_colour_3_tripod = MAXVAL;
#define RESETMARKS_VERTEX_COLOUR_3_TRIPOD {int mki; if ((markvalue_vertex_colour_3_tripod += 1) > MAXVAL) \
      { markvalue_vertex_colour_3_tripod = 1; for (mki=0;mki<MAX_VERTEX_COLOUR_DISTANCE_THREE_TRIPOD + 1;++mki) marks_vertex_colour_three_tripod[mki]=0;}}
#define MARK_VERTEX_COLOUR_3_TRIPOD(v) marks_vertex_colour_three_tripod[v] = markvalue_vertex_colour_3_tripod
#define ISMARKED_VERTEX_COLOUR_3_TRIPOD(v) (marks_vertex_colour_three_tripod[v] == markvalue_vertex_colour_3_tripod)

unsigned int marks_vertex_colour_four_tripod[MAX_VERTEX_COLOUR_DISTANCE_FOUR_TRIPOD + 1];
static int markvalue_vertex_colour_4_tripod = MAXVAL;
#define RESETMARKS_VERTEX_COLOUR_4_TRIPOD {int mki; if ((markvalue_vertex_colour_4_tripod += 1) > MAXVAL) \
      { markvalue_vertex_colour_4_tripod = 1; for (mki=0;mki<MAX_VERTEX_COLOUR_DISTANCE_FOUR_TRIPOD + 1;++mki) marks_vertex_colour_four_tripod[mki]=0;}}
#define MARK_VERTEX_COLOUR_4_TRIPOD(v) marks_vertex_colour_four_tripod[v] = markvalue_vertex_colour_4_tripod
#define ISMARKED_VERTEX_COLOUR_4_TRIPOD(v) (marks_vertex_colour_four_tripod[v] == markvalue_vertex_colour_4_tripod)


/* Some variables and marks used to test the colourability of graphs */
unsigned char number_of_colours_snarks[MAXN];
unsigned char colours_snarks[MAXN][REG];
unsigned char neighbour_index[MAXN][MAXN];

#define find_next_vertex(v, colour) (colours_snarks[v][0] == colour ? current_graph[v][0] : colours_snarks[v][1] == colour ? current_graph[v][1] : current_graph[v][2])

static int markvalue_snarks = MAXVAL;
unsigned int marks_snarks[MAXN][REG];
#define RESETMARKS_SNARKS {int mki, mkj; if ((markvalue_snarks += 1) > MAXVAL) \
      { markvalue_snarks = 1; for(mki=0;mki<MAXN;++mki) for(mkj=0;mkj<REG;++mkj) marks_snarks[mki][mkj]=0;}}
#define MARK_SNARKS(v, w) marks_snarks[v][w] = markvalue_snarks
#define UNMARK_SNARKS(v, w) marks_snarks[v][w] = markvalue_snarks - 1
#define ISMARKED_SNARKS(v, w) (marks_snarks[v][w] == markvalue_snarks)


/* Marks for the colour cycles */
static int markvalue_cycle_colour1 = MAXVAL;
unsigned int marks_cycle_colour1[MAXN][REG];
#define RESETMARKS_CYCLE_COLOUR1 {int mki, mkj; if ((markvalue_cycle_colour1 += 1) > MAXVAL) \
      { markvalue_cycle_colour1 = 1; for(mki=0;mki<MAXN;++mki) for(mkj=0;mkj<REG;++mkj) marks_cycle_colour1[mki][mkj]=0;}}
#define MARK_CYCLE_COLOUR1(v, w) marks_cycle_colour1[v][w] = markvalue_cycle_colour1
#define ISMARKED_CYCLE_COLOUR1(v, w) (marks_cycle_colour1[v][w] == markvalue_cycle_colour1)

static int markvalue_cycle_colour2 = MAXVAL;
unsigned int marks_cycle_colour2[MAXN][REG];
#define RESETMARKS_CYCLE_COLOUR2 {int mki, mkj; if ((markvalue_cycle_colour2 += 1) > MAXVAL) \
      { markvalue_cycle_colour2 = 1; for(mki=0;mki<MAXN;++mki) for(mkj=0;mkj<REG;++mkj) marks_cycle_colour2[mki][mkj]=0;}}
#define MARK_CYCLE_COLOUR2(v, w) marks_cycle_colour2[v][w] = markvalue_cycle_colour2
#define ISMARKED_CYCLE_COLOUR2(v, w) (marks_cycle_colour2[v][w] == markvalue_cycle_colour2)

static int markvalue_cycle_colour3 = MAXVAL;
unsigned int marks_cycle_colour3[MAXN][REG];
#define RESETMARKS_CYCLE_COLOUR3 {int mki, mkj; if ((markvalue_cycle_colour3 += 1) > MAXVAL) \
      { markvalue_cycle_colour3 = 1; for(mki=0;mki<MAXN;++mki) for(mkj=0;mkj<REG;++mkj) marks_cycle_colour3[mki][mkj]=0;}}
#define MARK_CYCLE_COLOUR3(v, w) marks_cycle_colour3[v][w] = markvalue_cycle_colour3
#define ISMARKED_CYCLE_COLOUR3(v, w) (marks_cycle_colour3[v][w] == markvalue_cycle_colour3)


unsigned long long int forbidden_edges_col_1_2_a[3*MAXN/2][3*MAXN/2];

static int markvalue_snarks_col_1_2_a = MAXVAL;
unsigned int marks_snarks_col_1_2_a[3*MAXN/2][3*MAXN/2];
#define RESETMARKS_SNARKS_COL_1_2_A {int mki, mkj; if ((markvalue_snarks_col_1_2_a += 1) > MAXVAL) \
      { markvalue_snarks_col_1_2_a = 1; for(mki=0;mki<3*MAXN/2;++mki) for(mkj=0;mkj<3*MAXN/2;++mkj) marks_snarks_col_1_2_a[mki][mkj]=0;}}
#define MARK_SNARKS_COL_1_2_A(v, w) marks_snarks_col_1_2_a[v][w] = markvalue_snarks_col_1_2_a
#define ISMARKED_SNARKS_COL_1_2_A(v, w) (marks_snarks_col_1_2_a[v][w] == markvalue_snarks_col_1_2_a)


unsigned long long int forbidden_edges_col_2_3_a[3*MAXN/2][3*MAXN/2];

static int markvalue_snarks_col_2_3_a = MAXVAL;
unsigned int marks_snarks_col_2_3_a[3*MAXN/2][3*MAXN/2];
#define RESETMARKS_SNARKS_COL_2_3_A {int mki, mkj; if ((markvalue_snarks_col_2_3_a += 1) > MAXVAL) \
      { markvalue_snarks_col_2_3_a = 1; for(mki=0;mki<3*MAXN/2;++mki) for(mkj=0;mkj<3*MAXN/2;++mkj) marks_snarks_col_2_3_a[mki][mkj]=0;}}
#define MARK_SNARKS_COL_2_3_A(v, w) marks_snarks_col_2_3_a[v][w] = markvalue_snarks_col_2_3_a
#define ISMARKED_SNARKS_COL_2_3_A(v, w) (marks_snarks_col_2_3_a[v][w] == markvalue_snarks_col_2_3_a)

unsigned long long int forbidden_edges_col_1_3_a[3*MAXN/2][3*MAXN/2];

static int markvalue_snarks_col_1_3_a = MAXVAL;
unsigned int marks_snarks_col_1_3_a[3*MAXN/2][3*MAXN/2];
#define RESETMARKS_SNARKS_COL_1_3_A {int mki, mkj; if ((markvalue_snarks_col_1_3_a += 1) > MAXVAL) \
      { markvalue_snarks_col_1_3_a = 1; for(mki=0;mki<3*MAXN/2;++mki) for(mkj=0;mkj<3*MAXN/2;++mkj) marks_snarks_col_1_3_a[mki][mkj]=0;}}
#define MARK_SNARKS_COL_1_3_A(v, w) marks_snarks_col_1_3_a[v][w] = markvalue_snarks_col_1_3_a
#define ISMARKED_SNARKS_COL_1_3_A(v, w) (marks_snarks_col_1_3_a[v][w] == markvalue_snarks_col_1_3_a)


unsigned long long int forbidden_edges_col_1_2_b[3*MAXN/2][3*MAXN/2];

static int markvalue_snarks_col_1_2_b = MAXVAL;
unsigned int marks_snarks_col_1_2_b[3*MAXN/2][3*MAXN/2];
#define RESETMARKS_SNARKS_COL_1_2_B {int mki, mkj; if ((markvalue_snarks_col_1_2_b += 1) > MAXVAL) \
      { markvalue_snarks_col_1_2_b = 1; for(mki=0;mki<3*MAXN/2;++mki) for(mkj=0;mkj<3*MAXN/2;++mkj) marks_snarks_col_1_2_b[mki][mkj]=0;}}
#define MARK_SNARKS_COL_1_2_B(v, w) marks_snarks_col_1_2_b[v][w] = markvalue_snarks_col_1_2_b
#define ISMARKED_SNARKS_COL_1_2_B(v, w) (marks_snarks_col_1_2_b[v][w] == markvalue_snarks_col_1_2_b)


unsigned long long int forbidden_edges_col_2_3_b[3*MAXN/2][3*MAXN/2];

static int markvalue_snarks_col_2_3_b = MAXVAL;
unsigned int marks_snarks_col_2_3_b[3*MAXN/2][3*MAXN/2];
#define RESETMARKS_SNARKS_COL_2_3_B {int mki, mkj; if ((markvalue_snarks_col_2_3_b += 1) > MAXVAL) \
      { markvalue_snarks_col_2_3_b = 1; for(mki=0;mki<3*MAXN/2;++mki) for(mkj=0;mkj<3*MAXN/2;++mkj) marks_snarks_col_2_3_b[mki][mkj]=0;}}
#define MARK_SNARKS_COL_2_3_B(v, w) marks_snarks_col_2_3_b[v][w] = markvalue_snarks_col_2_3_b
#define ISMARKED_SNARKS_COL_2_3_B(v, w) (marks_snarks_col_2_3_b[v][w] == markvalue_snarks_col_2_3_b)

unsigned long long int forbidden_edges_col_1_3_b[3*MAXN/2][3*MAXN/2];

static int markvalue_snarks_col_1_3_b = MAXVAL;
unsigned int marks_snarks_col_1_3_b[3*MAXN/2][3*MAXN/2];
#define RESETMARKS_SNARKS_COL_1_3_B {int mki, mkj; if ((markvalue_snarks_col_1_3_b += 1) > MAXVAL) \
      { markvalue_snarks_col_1_3_b = 1; for(mki=0;mki<3*MAXN/2;++mki) for(mkj=0;mkj<3*MAXN/2;++mkj) marks_snarks_col_1_3_b[mki][mkj]=0;}}
#define MARK_SNARKS_COL_1_3_B(v, w) marks_snarks_col_1_3_b[v][w] = markvalue_snarks_col_1_3_b
#define ISMARKED_SNARKS_COL_1_3_B(v, w) (marks_snarks_col_1_3_b[v][w] == markvalue_snarks_col_1_3_b)



/* Variables for nauty */
int lab[MAXN], ptn[MAXN], orbits[MAXN];
static DEFAULTOPTIONS_SPARSEGRAPH(options);
statsblk stats;
setword workspace[WORKSIZE];

sparsegraph sg; /* Sparse graph datastructure for nauty */
sparsegraph sg_canon; /* Sparse graph datastructure for nauty */

int use_defaultptn = 1;

int generators[MAXN][MAXN];
int number_of_generators;


/* Some auxiliary variables */
int max_edgepairlist_size[MAXN];
int max_edgepairlist_size_triangle = 0;

#define DEFAULT_MAX_EDGE4TUPLELIST_SIZE 1000
int max_edge4tuplelist_size = DEFAULT_MAX_EDGE4TUPLELIST_SIZE;

EDGE4TUPLE *edge_4tuples_list;

#define DEFAULT_MAX_EDGE4TUPLELIST_SIZE_G7 1000
int max_edge4tuplelist_size_g7 = DEFAULT_MAX_EDGE4TUPLELIST_SIZE_G7;

EDGE4TUPLE *edge_4tuples_list_g7;


#define MAX_NUMBER_OF_CYCLES (2*MAXN)
int current_cycle_number = 0;

//A bitvector of all cycles where the edge with label i is part of
setword edge_in_cycles[3 * MAXN / 2];

//The edgelist has to contain more than this number of elements, otherwise
//the colour cycles won't be searched
int min_edgepairlist_size = 0;


int binom_coefficients[MAXN][MAXN];

//Statistics
unsigned long long int nauty_calls = 0;

unsigned char *graphlist[MAXN + 1]; /* The list/buffer of graphs that still have to be written to the output file */
int graphlist_number_of_graphs[MAXN + 1]; /* The current number of graphs (of a certain order) in the graphlist */
int codelength[MAXN + 1]; /* The codelength of a graph of a certain order */
int listlength = 0; /* The maximum number of graphs that can be written in the buffer */


#define DEFAULT_MAX_GRAPHLIST_SNARKS_SIZE 100
//A list of snarkchildren of a certain parent
GRAPH *graphlist_snarks;
int graphlist_snarks_size = 0;
int max_graphlist_snarks_size = DEFAULT_MAX_GRAPHLIST_SNARKS_SIZE;

//The canonical graphs corresponding with the snarkchildren of graphlist_snarks
sparsegraph *canon_graphs;

//Sometimes delaying canonicity check in case of snarks
int still_has_to_check_if_graph_is_canonical = 0;

//Filename of the outputfile
char outputfilename[MAXN + 1][50];
FILE *outputfile[MAXN + 1];

//Filename of the outputfile for the irreducible graphs
char outputfilename_irred[50];
FILE *outputfile_irred;


/* Parameters of the program */
int singleout = 0; /* If singleout == 0, the program also writes graphs with <= number_of_vertices vertices */
int output_to_file = 1; /* If output_to_file == 1, the generated graphs are written to the file with name outputfilename, else they are written to stdout */
int graph6_output = 0; /* If graph6_output == 1, it will output the graphs in graph6 format, else it will output them in multicode format */
int noout = 0; /* If noout == 1, the generated graphs arent output, but only counted */
int output_prime_graphs = 0; /* If output_irred_graphs == 1, also the irreducible graphs will be written */

int check_cyclic_connectivity = 0;
int min_cyclic_connectivity = 0; /* The minimal cyclic edge connectivity */

int girth = 0; /* The minimum girth of the graphs which are generated */
int snarks = 0;

int modulo = 0;
int rest;
int mod;
int splitlevel; /* For modolo, determines at what level/depth of the recursion tree the calculation should be split */
unsigned int splitlevel_counter = 0;


/* Variables for girth == 6*/

//All pairs of disjoint pentagons!
#define MAX_NUM_PENTAGON_PAIRS (MAXN*MAXN) //Temp bound!!!
static unsigned long long int stored_disjoint_pentagon_pairs_edges[MAX_NUM_PENTAGON_PAIRS];
static int num_stored_disjoint_pentagon_pairs_edges = 0;

//All triples of disjoint pentagons!
#define MAX_NUM_PENTAGON_TRIPLES (MAXN*MAXN) //Temp bound!!!
static unsigned long long int stored_disjoint_pentagon_triples_edges[MAX_NUM_PENTAGON_TRIPLES];
static int num_stored_disjoint_pentagon_triples_edges = 0;

//Edge bitvectors of the pentagons where the inserted edge will be part of
#define MAX_NUM_FORBIDDEN_EDGES MAXN
unsigned long long int stored_forbidden_pentagon_edges[MAX_NUM_FORBIDDEN_EDGES];
int num_stored_forbidden_edges = 0;


//All pairs of hexagons which share at least 1 edge
#define MAX_NUM_HEXAGON_PAIRS_COMMON (MAXN*MAXN) //Temp bound!!!
static unsigned long long int stored_hexagon_pairs_edges_common[MAX_NUM_HEXAGON_PAIRS_COMMON];
static int num_stored_hexagon_pairs_edges_common = 0;

//All triples of hexagons which share at least 1 edge
#define MAX_NUM_HEXAGON_TRIPLES_COMMON (MAXN*MAXN) //Temp bound!!!
static unsigned long long int stored_hexagon_triples_edges_common[MAX_NUM_HEXAGON_TRIPLES_COMMON];
static int num_stored_hexagon_triples_edges_common = 0;

//All triples of hexagons which share at least 1 edge
#define MAX_NUM_HEXAGONS_PLUS_PENTAGON_TRIPLES_COMMON (MAXN*MAXN) //Temp bound!!!
static unsigned long long int stored_hexagons_plus_pentagon_triples_edges_common_hexagons[MAX_NUM_HEXAGONS_PLUS_PENTAGON_TRIPLES_COMMON];
static unsigned long long int stored_hexagons_plus_pentagon_triples_edges_common_pentagon[MAX_NUM_HEXAGONS_PLUS_PENTAGON_TRIPLES_COMMON];
static int num_stored_hexagons_plus_pentagon_triples_edges_common = 0;


//All pairs of heptagons which share at least 1 edge
#define MAX_NUM_HEPTAGON_PAIRS_COMMON (MAXN*MAXN) //Temp bound!!!
static unsigned long long int stored_heptagon_pairs_edges_common[MAX_NUM_HEPTAGON_PAIRS_COMMON];
static int num_stored_heptagon_pairs_edges_common = 0;

//All triples of heptagons which share at least 1 edge
#define MAX_NUM_HEPTAGON_TRIPLES_COMMON (MAXN*MAXN) //Temp bound!!!
static unsigned long long int stored_heptagon_triples_edges_common[MAX_NUM_HEPTAGON_TRIPLES_COMMON];
static int num_stored_heptagon_triples_edges_common = 0;

//All quadruples of heptagons which share at least 1 edge
#define MAX_NUM_HEPTAGON_QUADRUPLES_COMMON (2*MAXN*MAXN) //Temp bound!!!
static unsigned long long int stored_heptagon_quadruples_edges_common[MAX_NUM_HEPTAGON_QUADRUPLES_COMMON];
static int num_stored_heptagon_quadruples_edges_common = 0;


static int current_path[MAXN];
static int current_pathsize = 0;
static setword current_path_bitvector = 0;
static int contains_more_than_2_pentagons = 0;
static int contains_more_than_3_pentagons = 0;
static int contains_more_than_4_pentagons = 0;

int search_all_pentagons = 0;

//Temp bound!
#define MAX_NUM_PENTAGONS (10*MAXN)

static int disjoint_pentagon_index0;
static int disjoint_pentagon_index1;
static int disjoint_pentagon_index2;
static int disjoint_pentagon_index3;
static setword stored_pentagons[MAX_NUM_PENTAGONS];
static int num_pentagons_stored = 0;

static int stored_pentagons_vertices[MAX_NUM_PENTAGONS][MAXN];

//Warning: is limited to 42 vertices!
//Important: use lli instead of setword here so we can go up to 64 elements
static unsigned long long int stored_pentagons_edges[MAX_NUM_PENTAGONS];


static int contains_more_than_2_squares = 0;


//Temp bound!
#define MAX_NUM_SQUARES (10*MAXN)

static setword stored_squares[MAX_NUM_SQUARES];
static int num_squares_stored = 0;
static int stored_squares_vertices[MAX_NUM_SQUARES][MAXN];

//Warning: is limited to 42 vertices!
//Important: use lli instead of setword here so we can go up to 64 elements
static unsigned long long int stored_squares_edges[MAX_NUM_SQUARES];



//Temp bound!
#define MAX_NUM_HEXAGONS (10*MAXN)

static setword stored_hexagons[MAX_NUM_HEXAGONS];
static int num_hexagons_stored = 0;
//Warning: is limited to 42 vertices!
//Important: use lli instead of setword here so we can go up to 64 elements
static unsigned long long int stored_hexagon_edges[MAX_NUM_HEXAGONS];

static int stored_hexagon_vertices[MAX_NUM_HEXAGONS][MAXN];

//Temp bound!
#define MAX_NUM_HEPTAGONS (10*MAXN)

static setword stored_heptagons[MAX_NUM_HEPTAGONS];
static int num_heptagons_stored = 0;
//Warning: is limited to 42 vertices!
//Important: use lli instead of setword here so we can go up to 64 elements
static unsigned long long int stored_heptagon_edges[MAX_NUM_HEPTAGONS];


static setword hexagon_vertices = 0;

unsigned char nonhexagon_vertices[MAXN];
int num_nonhexagon_vertices = 0;

unsigned char nonhexagon_vertices_tripod[MAXN];
unsigned char nonhexagon_vertices_tripod_nbrs[MAXN][REG];
unsigned char nonhexagon_vertices_tripod_nbrs_size[MAXN];
int num_nonhexagon_vertices_tripod = 0;


//Only compute distance <= 3
#define MAX_DISTANCE 3

//Only compute distance <= 2
//#define MAX_DISTANCE_TRIPOD 2 //in case of girth 6
//#define MAX_DISTANCE_TRIPOD 1 //in case of girth 5

int max_distance_tripod;

//The distance from vertex i to vertex j
//Will only be computed for nonhexagon vertices
int distance[MAXN][MAXN];


int apply_tripod_optimisation = 0;
int real_girth = 0;
int test_for_snarks_tripod = 0;

int search_for_graphs_with_girth7 = 0;

int test_for_snarks_h_operation_girth7 = 0;

/********************************Methods***************************************/

/* Methods for the generation of irreducible graphs */
void generate_irreducible_graphs(int order, int min_girth, void (*userproc) (unsigned char (*)[REG + 1], int));
void init_generate_irreducible_graphs();

void extend_irreducible_graph(int edge_inserted);

void generate_eligible_diamond_edges(EDGE eligible_diamond_edges[], int *eligible_diamond_edges_size);
void edge_diamond_extend(EDGE eligible_diamond_edges[], int eligible_diamond_edges_size);

void generate_eligible_lollipop_edges(EDGE eligible_lollipop_edges[], int *eligible_lollipop_edges_size);
void edge_lollipop_extend(EDGE eligible_lollipop_edges[], int eligible_lollipop_edges_size);

void add_eligible_edge(unsigned char from, unsigned char to);
void replace_eligible_edge(int old_from, int old_to, int from, int to);

void add_nonadj_edge_diamond_to_list(unsigned char v0, unsigned char v1, unsigned char v2, unsigned char v3);

unsigned char determine_external_diamond_neighbour(IRRED_TRIANGLE diamond);
unsigned char determine_external_diamond_neighbour_index(IRRED_TRIANGLE diamond, int index);

void generate_non_adjacent_diamond_edge_pairs(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size);
void generate_nonadj_diamond_edge_pairs_one_lollipop(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size);
void generate_all_nonadj_diamond_edge_pairs(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size);

void nonadj_edge_diamond_extend(EDGEPAIR edge_pairs_list[], int edge_pair_list_size);

void update_edge_diamonds();

void update_irreducible_triangles_bitvector();
void fill_list_of_irreducible_triangles();


/* Methods for the generation of graphs */
void replace_neighbour(unsigned char vertex, unsigned char v1, unsigned char v2);

void extend(int edge_inserted, int trivial_group);
void edge_extend(EDGEPAIR edge_pairs_list[], int edge_pair_list_size);
void triangle_extend_all(int generators_local[][MAXN], int number_of_generators_local,
        int vertex_orbits_local[], int number_of_vertex_orbits, int groupsize);
void tripod_extend(EDGETRIPLE edge_triples_list[], int edge_triples_list_size);
void fourtuple_extend(EDGE4TUPLE edge_4tuples_list[], int edge_4tuples_list_size);

int is_part_of_irreducible_triangle_bitvector(int vertex);
int is_part_of_irreducible_triangle_diamond(int vertex, int *diamond);
int is_part_of_reducible_triangle(int vertex, int *triangle);

int is_a_bridge(unsigned char from, unsigned char to);
int is_a_bridge_list(unsigned char from, unsigned char to);


int contains_squares();
int find_squares(SQUARE squares[], setword squares_bitvectors[], int *squares_size, EDGE adjacent_squares[], int *adjacent_squares_sizes);
int can_destroy_pentagons_girth_at_least_4();
int are_adjacent_triangles();

int contains_more_than_3_disjoint_pentagons();
int contains_more_than_3_disjoint_squares();
static int contains_more_than_three_disjoint_pentagons();

void search_for_all_pentagons();

int edge_is_part_of_square(EDGE edge);
int inserted_edge_will_be_part_of_square(EDGEPAIR edgepair);

void find_cutvertices();
void determine_and_store_non_hex_neighbours(unsigned char vertex, int min_num_nonhex_neighbours);
void determine_and_store_non_pent_neighbours(unsigned char vertex, int min_num_nonhex_neighbours);

/* Methods to find the vertexsets */
void find_vertexsets(int vertices, unsigned char vertexset[][vertices], int *vertexset_size);
void generate_triangle_vertexsets_girth4(int num_vertices_in_set, unsigned char vertexset[][num_vertices_in_set], int *vertexset_size);
void find_vertexsets_girth4(int num_vertices_in_set, unsigned char vertexset[][num_vertices_in_set], int *vertexset_size);
void generate_all_vertexsets_girth4(int num_vertices_in_set, unsigned char vertexset[][num_vertices_in_set], int *vertexset_size);

void generate_all_vertexsets_triangles(int current_index, int num_vertices_in_set,
        unsigned char vertexset[][num_vertices_in_set], int *vertexset_size, unsigned char vertexset_temp[][num_vertices_in_set], int *vertexset_size_temp);
void generate_triangle_vertexsets(int current_index, int num_vertices_in_set, unsigned char vertexset[][num_vertices_in_set], int *vertexset_size);
void generate_triangle_vertexsets_remaining(int current_index, int num_vertices_in_set, unsigned char vertexset[][num_vertices_in_set], int *vertexset_size);
void generate_all_vertexsets(int current_index, int num_vertices_in_set, unsigned char vertexset[][num_vertices_in_set], int *vertexset_size);

int can_be_added_to_reducible_triangle_vertexset(unsigned char current_vertexset[], int current_size, int vertex, int triangle);

int is_part_of_same_reducible_triangle(int vertex, int triangle);
int is_part_of_same_irreducible_triangle(unsigned char vertex, int diamond);


/* Methods to find the edgepairs */
void find_edge_pairs(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size);
void generate_edgepairs_one_triangle(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size);
void generate_edgepairs_triangle_free_one_diamond(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size);

void generate_edgepairs_two_triangles(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size);

int generate_edgepairs_penultimate_level_girth4(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size);
int generate_edgepairs_penultimate_level_girth5_no_triangles(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size);
int generate_edgepairs_penultimate_level_girth6_no_squares(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size);

int find_squares_contains_diamonds(SQUARE squares[], setword squares_bitvectors[], int *squares_size, EDGE adjacent_squares[], int *adjacent_squares_sizes);
void generate_non_adjacent_edge_pairs_square_free(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size, int check_remove);

void generate_non_adjacent_edge_pairs_square_free_girth6(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size, int check_remove);
void generate_non_adjacent_edge_pairs_square_free_girth6_two_disjoint_pentagons(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size, int check_remove);
void generate_non_adjacent_edge_pairs_square_free_girth6_at_least_one_pentagon(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size, int check_remove);

void generate_non_adjacent_edge_pairs_square_free_contains_diamonds(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size);
void generate_edgepairs_no_triangles(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size);
void generate_non_adjacent_edge_pairs_one_square_contains_diamonds(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size, SQUARE square);

/* Methods for girth 5 */
void generate_non_adjacent_edge_pairs_one_square(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size, SQUARE square);
void generate_non_adjacent_edge_pairs_two_squares(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size, SQUARE square1, SQUARE square2);
void generate_non_adjacent_edge_pairs_two_adjacent_squares(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size, SQUARE squares[], EDGE common_edge);
void generate_non_adjacent_edge_pairs_four_squares(EDGEPAIR edge_pairs_list[], int *edge_pair_list_size, setword squares_bitvectors[], EDGE adjacent_squares[], int adjacent_squares_sizes);
void determine_common_edge(setword bitvector, EDGE common_edge);

int squares_have_multiple_common_neighbours(SQUARE squares[2]);


void transform_edge_into_canonical_form(EDGE edge);
void transform_edgepair_into_canonical_form(EDGEPAIR edgepair);
void transform_triangle_into_canonical_form(TRIANGLE triangle);
void transform_triangle_into_canonical_form_full(TRIANGLE triangle);
void transform_vertexset_into_canonical_form(unsigned char vertexset[], int size);
void transform_edgetriple_into_canonical_form(EDGETRIPLE edgetriple);
void transform_edgetriple_into_canonical_form_edges_and_first_pair_ok(EDGETRIPLE edgetriple);

void transform_vertexset_into_triangles(unsigned char vertexset[], int vertexset_size);
void undo_vertexset_triangle_transformation(unsigned char vertexset[], int vertexset_size);
int update_bridges_triangle_insert(unsigned char org_vertex, unsigned char new_vertex, unsigned char fixed_vertex);
int other_element_of_vertexset_is_part_of_diamond(IRRED_TRIANGLE diamond, unsigned char vertexset[], int vertexset_size, int index);

void remove_irreducible_triangle(int index);

void add_edge(EDGEPAIR edge_pair);
void remove_edge(EDGEPAIR edge_pair);
int contains_bridge(EDGEPAIR edge_pair);

void update_bridges_add_edge();
void add_bridge(unsigned char from, unsigned char to);
void remove_bridge(int index);
void replace_bridge(int old_from, int old_to, int from, int to);

void add_tripod(EDGETRIPLE edge_triple);
void remove_tripod(EDGETRIPLE edge_triple);

void add_4_tuple(EDGETRIPLE edge_4tuple);
void remove_4_tuple(EDGE4TUPLE edge_4tuple);

//Methods to calculate nauty partitions
void determine_major_edge_partitions();
void determine_vertex_partitions();

setword determine_vertex_neighbours_distance_three(unsigned char vertex);
setword determine_vertex_neighbours_distance_three_check_mark3(unsigned char vertex);
setword determine_vertex_neighbours_distance_two(unsigned char vertex);
setword determine_vertex_neighbours_distance_four(unsigned char vertex);

static int contains_neighbour_with_larger_label(unsigned char vertex, unsigned char candidate_neighbour);

int vertex_is_part_of_hexagon_blank(unsigned char vertex);

int is_only_vertex_from_irred_triangle(unsigned char vertexset[], int num_vertices_in_set, int vertexset_index, int *irred_triangle);
void determine_nauty_partitions_triangle(unsigned char vertexset[], int vertexset_size, int vertex_orbits_local[], int number_of_vertex_orbits);

/* Methods for the union find algorithm */
void determine_edgepair_orbits(EDGEPAIR edge_pairs_list[], int edge_pair_list_size, int edgepair_orbits[], int *number_of_orbits);
void determine_edgetriple_orbits(EDGETRIPLE edge_triples_list[], int edge_triples_list_size, int edgetriple_orbits[], int *number_of_orbits);
void determine_edge4tuple_orbits(EDGE4TUPLE edge_4tuples_list[], int edge_4tuples_list_size, int edge4tuple_orbits[], int *number_of_orbits);
void determine_edge_orbits(EDGE edgelist[], int edgelist_size, int edge_orbits[], int *number_of_orbits);
void determine_vertexset_orbits(int generators_local[][MAXN], int number_of_generators_local,
        int num_vertices_in_set, unsigned char vertexset[][num_vertices_in_set], int vertexset_size, int vertexset_orbits[], int *number_of_orbits);
void union_elements(int edgepair_orbits[], int root_orbits_size[], int *number_of_orbits, int a, int b);
int find_root(int edgepair_orbits[], int i);

int apply_special_generator(unsigned char org_vertexset[], int num_vertices_in_set, int diamond, unsigned char permutated_vertexset[]);
int is_only_vertex_from_diamond(unsigned char vertexset[], int num_vertices_in_set, int irred_triangle, int index);


/* Methods for (nauty) colours / canonicity */
int is_major_edge(int partitions_already_set, sparsegraph *sparsegraph_canon);
void determine_last_edge(sparsegraph sparse_graph_canon, int lab[], EDGE lastedge);

int is_major_tripod(int partitions_already_set, sparsegraph *sparsegraph_canon);

setword get_neighbours_distance_one(EDGE edge);
int new_edge_has_min_colour(EDGEPAIR edge_pair); //Lookahead
int new_edge_has_min_colour_no_squares(EDGEPAIR edge_pair); //Lookahead

int is_reducible_edge(EDGE edge);
int determine_edge_colour(EDGE edge);
int has_min_colour(EDGE inserted_edge, int *edge_inserted, EDGE neighbours[4]);
int has_min_colour_cycle(EDGE inserted_edge, int *edge_inserted);
int has_min_colour_cycle_girth_5(EDGE inserted_edge, int *edge_inserted);

int has_min_colour_tripod_girth5(int *edge_inserted);
int has_min_colour_tripod_girth6(int *edge_inserted);
int has_min_colour_tripod_girth7(int *edge_inserted);

int determine_edge_colour_expensive(EDGE edge);
int has_min_colour_distance_three(EDGE inserted_edge, int *edge_inserted, EDGE min_colour_edges[], int number_of_min_colour_edges);
int has_min_colour_distance_four(EDGE inserted_edge, int *edge_inserted, EDGE min_colour_edges[], int number_of_min_colour_edges);
int has_min_colour_combination(EDGE inserted_edge, int *edge_inserted, EDGE min_colour_edges[], int number_of_min_colour_edges);


/* Methods for snarks */
void determine_even_cycle(EDGE start_edge, EDGE colours, unsigned char cycle[], int *cycle_size);
void init_search_cycles();
void search_cycles();
void search_cycles_square(SQUARE squares[], int number_of_squares);

void init_is_colourable(unsigned char number_of_colours[]);
int is_colourable();
int modify_existing_colouring();
int is_colourable_other_colouring();
int is_colourable_other_colouring_startvertex(int startvertex);
int is_conflicting_colouring(unsigned char colours[][REG], int current_vertex, int colour);
void determine_available_colours(int used_colour, EDGE available_colours);
void determine_uncoloured_vertex(int vertex, int *uncoloured_vertex, int *missing_colour);
void unmark_colours(EDGE nonfree_labelled[], int nonfree_labelled_size);
int label_nonfree_choices(int current_vertex, EDGE nonfree_labelled[], int *nonfree_labelled_size);
int colour_next_free_choice(int number_of_coloured_edges);


/* Methods for girth == 6 */
static void determine_all_disjoint_pentagon_pairs();
int can_destroy_pentagons_girth_at_least_4_mark2();

int search_all_hexagons();
int contains_more_than_3_disjoint_hexagons();

static int contains_more_than_two_disjoint_pentagons();
int contains_disjoint_pentagons();

int are_two_cutvertices_method(unsigned char v1, unsigned char v2);

int is_major_edge_h_operation(int partitions_already_set, sparsegraph *sparsegraph_canon);
int contains_more_than_4_disjoint_hexagons();
int has_twocut(GRAPH graph);
void determine_last_edge_h_operation(sparsegraph sparse_graph_canon, int lab[], EDGE lastedge);
int is_reducible_edge_H_operation(EDGE edge);
int has_min_colour_H_operation_girth7(EDGE inserted_edge, int *edge_inserted);
int has_min_colour_H_operation_girth6(EDGE inserted_edge, int *edge_inserted);
int colour_next_free_choice_from_zero(int number_of_coloured_edges);
int contains_more_than_2_disjoint_squares();

/* Auxiliary methods */
void code_multicode(unsigned char *code, GRAPH g, int num_of_vertices);

void aufschreiben();
void aufschreiben_already_filtered(GRAPH g, int num_vertices);
void wegspeichern(unsigned char *liste[MAXN], int num_of_vertices);

void aufschreiben_irred();
void wegspeichern_irred(unsigned char *liste, int codelength, int num_of_vertices);

void calculate_binom_coefficients(int max_n);
int power(int base, int exponent);

void printgraph_nauty(graph *g, int current_number_of_vertices);
void print_sparse_graph_nauty(sparsegraph sparse_graph);
void printgraph();
void printgraph_irred();

void init_nauty_options();
void save_generators(int count, int perm[], int orbits[],
        int numorbits, int stabvertex, int n);
void copy_sparse_graph();

/**
 * User defined method which is called each time a graph is generated
 * The result is a REG-regular graph. The neighbour lists only contain REG actual elements,
 * but using REG + 1 since it is slightly more efficient
 */
//void handle_3_regular_result(unsigned char snarkhunter_graph[MAXN][REG + 1], int order);
void (*userproc_sh) (unsigned char (*)[REG + 1], int) = NULL;


#endif	/* _SNARKHUNTER_H */
