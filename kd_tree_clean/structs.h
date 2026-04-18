#ifndef STRUCTS_H
#define STRUCTS_H

#include <omp.h>
#include <vector>
#include <semaphore.h>
//Structs


//kd tree node for example program
struct NodeEx{
    double * point;
    //used to keep track of which dimension the node should be sorted in based on the cd
    unsigned int pid; //id of the point
    struct NodeEx * left;
    struct NodeEx * right;
    unsigned int cnt;
};

//kd-tree node
struct Node{
    DTYPE * point;
    //used to keep track of which dimension the node should be sorted in based on the cd
    unsigned int pid; //id of the point
    struct Node * left;
    struct Node * right;
    unsigned int cnt;
};

//kd-tree node but with indices to children instead of pointers
struct NodeIndices{
    DTYPE * point;
    //used to keep track of which dimension the node should be sorted in based on the cd
    unsigned int pid; //id of the point
    unsigned int cnt;
    int leftIdx; //"NULL" pointer is -1 so need signed int
    int rightIdx; //"NULL" pointer is -1 so need signed int
};

//Each query rectangle is defined by the lower left and upper right
//points that define the minimum bounding box
//Also used for partitioning the space during range query searches
struct Rectangle{
    DTYPE min[IDIM];
    DTYPE max[IDIM];
    
};

struct RectangleEx{
    double min[IDIM];
    double max[IDIM];
    
};

//used for sorting the processing order of subdatasets
struct keyValPair{
    unsigned int key;
    unsigned int value;
};

//used for sorting the result set
struct keyValResultSet{
    unsigned int key;
    unsigned int value;
};

//used for sorting
struct keyVal{
DTYPE pointSortedDim;
unsigned int pid;
unsigned int index; //0, 1, 2,...,
};

struct CompareBySortedCoord {
    bool operator()(const keyVal& left, DTYPE right) const {
        return left.pointSortedDim < right;
    }

    bool operator()(DTYPE left, const keyVal& right) const {
        return left < right.pointSortedDim;
    }
};


//Used to sort the points based on the cutting dimension
struct keyValData{
DTYPE * point;
unsigned int pid;
unsigned int sortDim;
};

struct keyValueResultBuffers{
    unsigned int * keys;
    unsigned int * values;
};


struct resultSetStorage{
    //for COUNTONLY==0 and RESULTSETSTORAGE==1
    unsigned int ** resultKeys; //NUMSUBDATASETS number of pointers
    unsigned int ** resultValues; //NUMSUBDATASETS number of pointers
    unsigned int * resultSetSizeEachKeyValueArray; //size of NUMSUBDATASETS
    unsigned int NUMSUBDATASETS;
    //for COUNTONLY==0 and RESULTSETSTORAGE==2
    std::vector<unsigned int> *resultSetVect;
    //for COUNTONLY==1 (only count the number of points within epsilon)
    unsigned int * resultSetCounts;
};

// Used to store candidate set of points
// This is simply an array with a counter
struct Candidates{
    unsigned int candidates[CANDBUFFERSIZE];
    unsigned int cnt;
};


struct dim_reorder_sort
{
        unsigned int dim; //point dimension
        double variance; //variance of the points in this dimension
};


struct stats{
    double totalEndToEndTime;
    double totalIndexTime;
    // double totalDataTransferTime;
    // double totalDPUKernelOnlyComputeTime;
    // pthread_mutex_t statslock;
};




#endif