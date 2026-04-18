#include "params.h"

#include <unistd.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <cstdint>
#include <algorithm>
#include <cstring>
#include <cassert>
#include <numeric>
#include <limits.h>
#include <parallel/algorithm>
#include <algorithm>

#include "structs.h"
#include "compare.h"
#include "kdtree.h"



#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)


using namespace std;


//init struct members:
void initStatsStruct(stats * stats){
    stats->totalEndToEndTime=0.0;
    stats->totalIndexTime=0.0;
}


void printStatsKDTree(unsigned int NPOINTS, unsigned int NQPOINTS, stats * stats, char * inputDatasetFname, char * inputQuerysetFname, double epsilonDouble, uint64_t totalNeighbors)
{
    char fname[]="stats.txt";
    ofstream outputStats;
    outputStats.open(fname,ios::app); 

    double totalTimeWithoutIndex = (stats->totalEndToEndTime)-(stats->totalIndexTime);

    outputStats<<"Dataset/Queryset: "<< inputDatasetFname<<", "<< inputQuerysetFname<<", "<<"NPOINTS/NQPOINTS: "<<NPOINTS<<", "<<NQPOINTS<<", "<<epsilonDouble<<", "<<totalNeighbors<<", "<<\
    "Times (s): /Total/Total-Index/Index: "<<stats->totalEndToEndTime<<", "<<totalTimeWithoutIndex<<", "<<stats->totalIndexTime<<", "<<\
    "NUMCPUTHREADS/IDIM/DTYPE: "<<NUMCPUTHREADS<<", "<<IDIM<<", "<<STR(DTYPE)<<", "<<\
    "SHORTCIRCUIT: "<<SHORTCIRCUIT<<", "<<\
    "RESULTSETSTORAGECPU: "<<RESULTSETSTORAGECPU<<", "<<"COUNTONLY: "<<COUNTONLY<<endl;

    outputStats.close();
}





int importDataset(char * fname, unsigned int N, unsigned int NDIM, double * dataset)
{

    FILE *fp = fopen(fname, "r");

    if (!fp) {
        printf("Unable to open file\n");
        exit(0);
    }

    char buf[4096];
    int rowCnt = 0;
    int colCnt = 0;
    while (fgets(buf, 4096, fp) && rowCnt<N) {
        colCnt = 0;

        char *field = strtok(buf, ",");
        double tmp;
        sscanf(field,"%lf",&tmp);
        dataset[rowCnt*NDIM+colCnt]=tmp;

        
        while (field) {
          colCnt++;
          field = strtok(NULL, ",");
          
          if (field!=NULL)
          {
          double tmp;
          sscanf(field,"%lf",&tmp);
          dataset[rowCnt*NDIM+colCnt]=tmp;
          }   

        }
        rowCnt++;
    }

    fclose(fp);

    return 0;


}

keyValData * createKeyValStruct(DTYPE * dataset, unsigned int NPOINTS, unsigned int NDIM)
{
    keyValData * tmpDataset = (keyValData*)malloc(sizeof(keyValData)*NPOINTS);

    //populate struct
    for (unsigned int i=0; i<NPOINTS; i++){
    //cd dimension
    tmpDataset[i].point=&dataset[i*NDIM];
    //point id
    tmpDataset[i].pid=i;
    tmpDataset[i].sortDim=0;
    }


    return tmpDataset;
}



void convertDatasetDatatype(DTYPE * dataset, double * datasetDouble, unsigned int NPOINTS, unsigned int NDIM)
{
    for (uint64_t i=0; i<NPOINTS*NDIM; i++)
    {
        dataset[i]=(DTYPE)datasetDouble[i];
    }
}




//resultSetVectorsIn is of size NUMTHREADS
void populateResultVectorsForNeighbortableCPU(std::vector<unsigned int> * resultSetVectorsOut, 
    std::vector<unsigned int> * resultSetVectorsIn, const unsigned int NUMTHREADS)
{
    for(uint64_t i=0; i<NUMTHREADS; i++)
    {
        for(uint64_t j=0; j<(resultSetVectorsIn[i]).size(); j+=2)
        {
            unsigned int key = resultSetVectorsIn[i][j];
            unsigned int value = resultSetVectorsIn[i][j+1];
            resultSetVectorsOut[key].push_back(value);
        }
    }
}


void printNeighborTableCount(unsigned int * resultSetCounts, uint64_t pointsToPrint, uint64_t NPOINTS)
{
    for(uint64_t i=0; i<pointsToPrint; i++){
        printf("\n%lu, %u", i, resultSetCounts[i]);
    }
}

void printNeighborTable(std::vector<unsigned int> *resultSet, uint64_t pointsToPrint, uint64_t NPOINTS)
{

    uint64_t sum=0;
    for(uint64_t i=0; i<NPOINTS; i++){
        sum+=resultSet[i].size();
    }
    printf("\nSanity check: sum of result vector sizes (should be same as the total result set size): %lu", sum);

    printf("\nPrinting neighbortable (sorted): ");
    for(uint64_t i=0; i<pointsToPrint; i++){
        printf("\nPoint %lu: ", i);
        std::sort(resultSet[i].begin(), resultSet[i].end());
        for(unsigned int j=0; j<resultSet[i].size(); j++){
            printf("%u, ", resultSet[i][j]);
        }
    }
}

//if we reorder the dataset we also need to reorder the queryset
void ReorderByDimension(double * dataset, double * queryset, unsigned int NPOINTS, unsigned int NQPOINTS, unsigned int NDIM)
{

    double tstarttotal = omp_get_wtime();

    double sums[IDIM];
    double average[IDIM];
    struct dim_reorder_sort dim_variance[IDIM];
    for (int i=0; i< IDIM; i++){
        sums[i]=0;
        average[i]=0;
    }

    double greatest_variance=0;
    int greatest_variance_dim=0;

    
    int sample=100;
    double inv_sample=1.0/(sample*1.0);
    printf("\nCalculating variance based on on the following fraction of pts: %f",inv_sample);
    double tvariancestart=omp_get_wtime();
        //calculate the variance in each dimension  
        for (int i=0; i<IDIM; i++)
        {
            //first calculate the average in the dimension:
            //only use every 10th point
            for (int j=0; j<NPOINTS; j+=sample)
            {
            // sums[i]+=(*NDdataPoints)[j][i];
            sums[i]+=dataset[j*NDIM+i];
            }


            average[i]=(sums[i])/(NPOINTS*inv_sample);
            // printf("\nAverage in dim: %d, %f",i,average[i]);

            //Next calculate the std. deviation
            sums[i]=0; //reuse this for other sums
            for (int j=0; j<NPOINTS; j+=sample)
            {
            sums[i]+=(dataset[j*NDIM+i]-average[i])*(dataset[j*NDIM+i]-average[i]);
            }
            
            dim_variance[i].variance=sums[i]/(NPOINTS*inv_sample);
            dim_variance[i].dim=i;
            
            // printf("\IDIM:%d, variance: %f",dim_variance[i].dim,dim_variance[i].variance);

            if(greatest_variance<dim_variance[i].variance)
            {
                greatest_variance=dim_variance[i].variance;
                greatest_variance_dim=i;
            }
        }

    double tvarianceend = omp_get_wtime();
    

    
    double tqsortstart = omp_get_wtime();
    qsort(dim_variance, IDIM, sizeof(dim_reorder_sort), compareByDimVariance);    
    double tqsortend = omp_get_wtime();
    


    for (int i=0; i<IDIM; i++){
        printf("\nReodering dimension by: dim: %d, variance: %f",dim_variance[i].dim,dim_variance[i].variance);
    }

    printf("\nDimension with greatest variance: %d",greatest_variance_dim);

    //If there are only 2 dimensions and the dimension with the greatest variance is dimension 0, then we do not need to reorder by dimension
    //return and do not modify arrays
    if(IDIM==2 && greatest_variance_dim==0){
        printf("\nDimension 0 already has the greatest variance, disregarding reordering these dimensions");
        return;
    }

    //If there are only 2 dimensions and the first two dimensions have very similar variance (within 90%) but variance(d=0) < variance(d=1)
    //return and do not modify arrays
    if(IDIM==2 && greatest_variance_dim==1){ 
        double varianceRatio = dim_variance[1].variance/dim_variance[0].variance;
        if(varianceRatio>=0.90){
            printf("\nDimensions 0 and 1 have similar variance (ratio: %f), disregarding reordering these dimensions", varianceRatio);
            return;
        }
        
    }


    double tmemcpystart = omp_get_wtime();

    //copy the dataset
    double * tmpDataset = (double *)malloc(sizeof(double)*(uint64_t)NPOINTS*(uint64_t)NDIM);  
    //copy dataset into temp dataset
    // memcpy(tmpDataset, dataset, sizeof(double)*NPOINTS*NDIM);
    #pragma omp parallel for num_threads(8)
    for(uint64_t i=0; i<NPOINTS*NDIM; i++){
        tmpDataset[i] = dataset[i];
    }

    //copy the queryset
    double * tmpQueryset = (double *)malloc(sizeof(double)*(uint64_t)NQPOINTS*(uint64_t)NDIM);  
    //copy queryset into temp queryset
    // memcpy(tmpQueryset, queryset, sizeof(double)*NQPOINTS*NDIM);

    #pragma omp parallel for num_threads(8)
    for(uint64_t i=0; i<NQPOINTS*NDIM; i++){
        tmpQueryset[i] = queryset[i];
    }

    double tmemcpyend = omp_get_wtime();
    

    
    double treordercopystart = omp_get_wtime();
    
    
    //rewrite dimensions based on variance
    // for (unsigned int j=0; j<IDIM; j++){

    //     unsigned int originDim=dim_variance[j].dim;  
    //     //dataset:
    //     // #pragma omp parallel for num_threads(8) schedule(guided)
    //     for (unsigned int i=0; i<NPOINTS; i++)
    //     {   
    //         dataset[i*NDIM+j]=tmpDataset[i*NDIM+originDim];
    //     }

    //     //queryset:
    //     // #pragma omp parallel for num_threads(8) schedule(guided)
    //     for (unsigned int i=0; i<NQPOINTS; i++)
    //     {   
    //         queryset[i*NDIM+j]=tmpQueryset[i*NDIM+originDim];
    //     }
    // }

    #pragma omp parallel for num_threads(8)
    for (unsigned int i=0; i<NPOINTS; i++){
    
        // unsigned int originDim=dim_variance[j].dim;  
        //dataset:
    
        for (unsigned int j=0; j<IDIM; j++){     
            dataset[i*NDIM+j]=tmpDataset[i*NDIM+dim_variance[j].dim];
        }
    }


     //queryset:
        #pragma omp parallel for num_threads(8)
        for (unsigned int i=0; i<NQPOINTS; i++)
        {   
            // unsigned int originDim=dim_variance[j].dim;  
            for (unsigned int j=0; j<IDIM; j++){
            queryset[i*NDIM+j]=tmpQueryset[i*NDIM+dim_variance[j].dim];
            }
        }

    double treordercopyend = omp_get_wtime();

    double tendtotal = omp_get_wtime();

    double totalTime = tendtotal - tstarttotal;

    printf("\nTime to compute variance: %f (frac: %f)", (tvarianceend - tvariancestart), (tvarianceend - tvariancestart)/totalTime);
    printf("\nTime to memcpy: %f (frac: %f)", (tmemcpyend - tmemcpystart), (tmemcpyend - tmemcpystart)/totalTime);
    printf("\nTime to qsort: %f (frac: %f)", (tqsortend - tqsortstart), (tqsortend - tqsortstart)/totalTime);
    printf("\nTime to reorder/copy: %f (frac: %f)", (treordercopyend - treordercopystart), (treordercopyend - treordercopystart)/totalTime);
    printf("\nTotal Time to sort by variance: %f", totalTime);
    
}






int main(int argc, char *argv[])
{

    //Read in parameters from file:
    if (argc<7)
    {
    printf("\nIncorrect number of input parameters.  \nShould be dataset file, queryset file, number of datapoints, number of querypoints, epsilon, number of dimensions\n\n");
    return 0;
    }

    ///////////////////////////////////////////
    //start of reading data and input parameters
    

    //set OpenMP threads
    omp_set_num_threads(NUMCPUTHREADS);
    omp_set_max_active_levels(2);


    
    //copy parameters from commandline:
    char inputDatasetFname[500];
    char inputQuerysetFname[500];
    unsigned int NPOINTS; //datapoints |D|
    unsigned int NQPOINTS; //query points |Q|
    double epsilonDouble;
    unsigned int NDIM;
    
    

    strcpy(inputDatasetFname,argv[1]);
    strcpy(inputQuerysetFname,argv[2]);
    sscanf(argv[3],"%u",&NPOINTS);
    sscanf(argv[4],"%u",&NQPOINTS);
    sscanf(argv[5],"%lf",&epsilonDouble);
    sscanf(argv[6],"%u",&NDIM);


    if(STR(DTYPE)!="float" && STR(DTYPE)!="double"){
        printf("\nERROR: DTYPE is not float or double\n\n");
        return 0;
    }



       
    printf("\nDataset: %s", inputDatasetFname);
    printf("\nQueryset: %s", inputQuerysetFname);
    printf("\nNumber of data points: %u", NPOINTS);
    printf("\nNumber of query points: %u", NQPOINTS);
    printf("\nEpsilon (as double): %lf", epsilonDouble);
    printf("\nNumber of data dimensions: %u", NDIM);
    
    
    



    DTYPE epsilon;




    //dataset
    double * datasetDouble = (double*)malloc(sizeof(double)*NPOINTS*NDIM);
    importDataset(inputDatasetFname, NPOINTS, NDIM, datasetDouble);
    
    //queryset
    double * querysetDouble = (double*)malloc(sizeof(double)*NQPOINTS*NDIM);
    importDataset(inputQuerysetFname, NQPOINTS, NDIM, querysetDouble);


    DTYPE * dataset = (DTYPE*)malloc(sizeof(DTYPE)*NPOINTS*NDIM);
    DTYPE * queryset = (DTYPE*)malloc(sizeof(DTYPE)*NQPOINTS*NDIM);


    
    double timeReorder = 0.0;
    #if REORDERDIMS==1
    double tstartReorderByDimension = omp_get_wtime();
    //ReorderByDimension uses doubles
    ReorderByDimension(datasetDouble, querysetDouble, NPOINTS, NQPOINTS, NDIM);
    double tendReorderByDimension = omp_get_wtime();
    timeReorder = tendReorderByDimension - tstartReorderByDimension;
    printf("\nTime to reorder by dimension: %f", timeReorder);
    #endif



    if(STR(DTYPE)!="double"){
        printf("\nDTYPE is not the default (double). DTYPE is %s.", STR(DTYPE));

            
        
            //convert double to float (straightforward cast each element)
            if(STR(DTYPE)=="float"){
            epsilon = (float)epsilonDouble;
            convertDatasetDatatype(dataset, datasetDouble, NPOINTS, NDIM);
            convertDatasetDatatype(queryset, querysetDouble, NQPOINTS, NDIM);
            }

            

    }
    //if the datatype is double
    //do stupid copy instead of updating pointers because
    //of compiler
    else{
        epsilon = epsilonDouble;

        memcpy(dataset, datasetDouble, sizeof(DTYPE)*NPOINTS*NDIM);
        memcpy(queryset, querysetDouble, sizeof(DTYPE)*NQPOINTS*NDIM);

    }


    //free the double precision input dataset/queryset
    free(datasetDouble);
    free(querysetDouble);
    


    

    //end of reading data and input parameters
    ///////////////////////////////////////////
    

    double tstartIndex, tendIndex, timeIndex;

    
    //Inserts all nodes at once with a single recursive call
    //The tree is balanced and requires balance because the inserts do not start at the root each time
    //The count inside each node will be incorrect because the insertions do not start at the root
    //Uses indices into array instead of pointers in the tree structure
    
    ///////////////////////////////////////////
    //Start of kd-tree code

        
    //Flag to allow duplicates in the tree
    bool allowDuplicates = true;

    

    //stats output to file, stored/updated in struct
    stats stats;
    initStatsStruct(&stats);

    //index into arrays instead of pointer chasing
    int root;

    
    tstartIndex = omp_get_wtime();

    unsigned int startingDepth=0;
    //These are statistics on tree depth so we can determine how balanced it is
    unsigned int insertedNodeDepth=0;
    unsigned int maxDepth=0;

    //Create nodes for each point in the dataset so
    //that we do not allocate one at a time
    struct NodeIndices * treeNodes = (NodeIndices*)malloc(sizeof(NodeIndices)*NPOINTS);
    initPointsInPreallocatedNodes(treeNodes, dataset, NPOINTS, NDIM);

    //Create a struct that represents the dataset as key/value pairs so that we can sort on the cutting
    //dimension and keep track of the point ids
    keyValData * datasetAsKeyVal = createKeyValStruct(dataset, NPOINTS, NDIM);

    //Inserts all of the data points, but cannot be used for counting range queries (which is fine)
    insertedNodeDepth = 0;
    

    root = insertAllWithIdx(root, treeNodes, datasetAsKeyVal, NPOINTS, startingDepth, &maxDepth, &insertedNodeDepth, NDIM, allowDuplicates);
    // maxDepth = max(insertedNodeDepth, maxDepth);

    free(datasetAsKeyVal);

    printf("\nRoot idx: %u", root);
    printf("\nNPOINTS: %u", NPOINTS);
    printf("\nOptimal tree depth ceil(log_2(|D|): %0.3f, Maximum tree depth: %u", ceil(log2((double)NPOINTS)), maxDepth);
    printf("\nFor a balanced kd-tree the optimal depth should be 1 level higher than the maximum tree depth because the latter counts the root as level 0 and not level 1.");


    tendIndex = omp_get_wtime();
    timeIndex = tendIndex - tstartIndex;

    //Search and refine
    double tstart=omp_get_wtime();

    ///////////////////////////////////////
    //result set options: 

    //MG: April 15
    //I left all of these result set options here as they can be slightly modified to work with 
    //arrays instead of vectors in case vectors aren't possible on the FPGA

    //We do not store the candidate set
    #if COUNTONLY==0 && RESULTSETSTORAGECPU==0
    std::vector<unsigned int> *resultSet = NULL;
    #endif

    //One vector per thread that searches and refines
    //it stores both keys and values in the same array
    //elem 0: key, elem 1: value, and so on
    #if COUNTONLY==0 && RESULTSETSTORAGECPU==1
    std::vector<unsigned int> *resultSet = new vector<unsigned int>[NUMCPUTHREADS]; 
    #endif

    //one vector per query point
    //the index of the result set refers to the query point index
    #if COUNTONLY==0 && RESULTSETSTORAGECPU==2
    std::vector<unsigned int> *resultSet = new vector<unsigned int>[NQPOINTS]; 
    #endif

    //if we only perform counting range queries
    #if COUNTONLY==1
    unsigned int *resultSet = (unsigned int*)malloc(sizeof(unsigned int)*NQPOINTS); 
    #endif
    ///////////////////////////////////////

    //1- refers to using the tree with nodes of type NodeIndices
    //need the dataset for the global bounding volume
    uint64_t totalNeighbors = searchAndRefineWithMode(1, NULL, root, treeNodes, dataset, queryset, NDIM, NPOINTS, NQPOINTS, epsilon, resultSet);
    double tend=omp_get_wtime();
    double timeSearchAndRefine = tend - tstart;




    stats.totalEndToEndTime=timeIndex+timeSearchAndRefine+timeReorder;
    stats.totalIndexTime=timeIndex;

    printf("\n*******************************\n");
    printf("\nTime to index: %f", stats.totalIndexTime);
    printf("\nTotal Time to Search and Refine: %f", timeSearchAndRefine);
    printf("\nTotal Time including index construction: %f", stats.totalEndToEndTime);


    
    

    printStatsKDTree(NPOINTS, NQPOINTS, &stats, inputDatasetFname, inputQuerysetFname, epsilonDouble, totalNeighbors);



    //MG: April 15, 2026
    //For printing the results using the several methods above
    //this may be important if you cannot store the results in a vector so I did not remove
    //all of these options
    //None of these will run because they are set to if(False) but if you want to see the results printed
    //to the screen, change these to if(true)
    

    #if COUNTONLY==0 && RESULTSETSTORAGECPU==2

    if(false)
    {        
    //Print neighbor table result set:
    printNeighborTable(resultSet, NQPOINTS, NQPOINTS);

    //only first the first 1000 points
    // printNeighborTable(resultSet, 1000, NQPOINTS);
    
    }

        
    #endif

    #if COUNTONLY==0 && RESULTSETSTORAGECPU==1
    
    if(false)
    {
        //The result set is stored as key/values in vectors, with one vector per CPU thread
        //so to construct the neighbortable we need to do a few extra steps compared to RESULTSETSTORAGE==2
        //(FOR VALIDATION ONLY)
        
        std::vector<unsigned int> * resultSetVectors = new vector<unsigned int>[NQPOINTS]; 
        populateResultVectorsForNeighbortableCPU(resultSetVectors, resultSet, NUMCPUTHREADS);


        
        //Print neighbor table result set:
        //Print all (default)
        printNeighborTable(resultSetVectors, NQPOINTS, NQPOINTS);
        
        //Print first 1000 points
        // printNeighborTable(resultSetVectors, 1000, NQPOINTS);

        delete [] resultSetVectors;
    }

    #endif

    //printing the result of the counting range queries
    #if COUNTONLY==1

    if(false)
    {
        //Print all (default)
        printNeighborTableCount(resultSet, NQPOINTS, NQPOINTS);

        //Print first 1000 points
        // printNeighborTableCount(resultSet, 1000, NQPOINTS);
    }
    
    #endif


    
    printf("\n\n");
    return 0;

}

