//0- similarity searches that generate result set pairs
//1- counting similarity searches that only return the number of neighbors for each query point
#define COUNTONLY 0

//DTYPE- float or double
#define DTYPE double

//Optimization: Short circuit the distance calculation (only useful for moderate and higher dimensional datasets i.e., when the dimensionality > 4 or so)
#define SHORTCIRCUIT 1

//Optimization: Reorder dataset by variance (only useful for higher dimensional dataset i.e., when the dimensionality > IDIM)
#define REORDERDIMS 1

//For searching the kd-tree and refining and index construction
#define NUMCPUTHREADS 32 //32

//The number of indexed dimensions
//Maximum should likely be 6 dimensions or so
#define IDIM 2

//Result set storage for the CPU:
//Three options (default 2)
#define RESULTSETSTORAGECPU 2 
//0- do not store results
//1- each CPU thread stores key/value pairs in their own array 
//2- store result pairs in vectors where each vector corresponds to one query point

//number of points stored for each candidate set
//Need one of these buffers per thread for the kd-tree implementation
#define CANDBUFFERSIZE 1000000