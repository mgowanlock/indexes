//k-d tree implementation based on Prof. David Mount's excellent lecture notes here:
//https://www.cs.umd.edu/class/fall2019/cmsc420-0201/Lects/lect13-kd-dict.pdf
//https://www.cs.umd.edu/class/fall2019/cmsc420-0201/Lects/lect14-kd-query.pdf



//left partition of a point using the cutting dimension
//splits the boundingVolume rect into a leftpart(ition) and rightpart(ition)
Rectangle leftpart(unsigned int cd, DTYPE * point, Rectangle * boundingVolume)
{

    Rectangle tmpRect;

    memcpy(tmpRect.min, boundingVolume->min, sizeof(DTYPE)*IDIM);
    memcpy(tmpRect.max, boundingVolume->max, sizeof(DTYPE)*IDIM);

    //update the maximum based on the cutting dimension of the point
    tmpRect.max[cd]=point[cd];

    return tmpRect;

}

//right partition of a point using the cutting dimension
//splits the boundingVolume rect into a leftpart(ition) and rightpart(ition)
Rectangle rightpart(unsigned int cd, DTYPE * point, Rectangle * boundingVolume)
{

    Rectangle tmpRect;

    memcpy(tmpRect.min, boundingVolume->min, sizeof(DTYPE)*IDIM);
    memcpy(tmpRect.max, boundingVolume->max, sizeof(DTYPE)*IDIM);

    //update the maximum based on the cutting dimension of the point
    tmpRect.min[cd]=point[cd];

    return tmpRect;

}



void initNode(NodeIndices * treeNode)
{

    //when not initializing the point data
    treeNode->leftIdx = -1;
    treeNode->rightIdx = -1;
    // treeNode->sortDim = 0;
    //initialize the number of points in the node to be equal to 1 (the one we just added)
    treeNode->cnt = 1; 
}


bool checkEqualPoints(DTYPE * pointA, DTYPE * pointB)
{
    for(unsigned int i=0; i<IDIM; i++){
        if(pointA[i]!=pointB[i]){
            return false;
        }
    }
    return true;
}

void generateQueryMBBFromPoint(Rectangle * queryMBB, DTYPE * point, DTYPE epsilon)
{
    for (unsigned int i=0; i<IDIM; i++){
    queryMBB->min[i]=point[i]-epsilon;
    queryMBB->max[i]=point[i]+epsilon;
    }
}




//returns true if the point is contained within the MBR
bool rangeQueryPointContained(Rectangle * rangeQueryMBR, Node * root)
{
    for (unsigned int i=0; i<IDIM; i++)
    {
        if(!(root->point[i]>=rangeQueryMBR->min[i] && root->point[i]<=rangeQueryMBR->max[i]))
        {
            return false;
        }

    }
    return true;
}

void updateCandidateSet(Node * root, Candidates * candSet)
{

    if(candSet->cnt==CANDBUFFERSIZE){
        printf("\n[Error] Candidate set overflow! Candidate set size is %u", candSet->cnt);
        exit(0);
    }

    //add buffer overflow logic
    //add current buffer size to the candSet struct as well 
    //(reinventing vectors)
    candSet->candidates[candSet->cnt] = root->pid;
    candSet->cnt++;
}


//Determines if the bounding volume is contained within the range query MBB
// bool rangeQueryBoundingVolumeContained(Rectangle * rangeQueryMBR, Rectangle * boundingVolume)
// {
    
//     for (unsigned int i=0; i<IDIM; i++)
//     {
//         if(!((rangeQueryMBR->min[i] <= boundingVolume->min[i]) && 
//             (rangeQueryMBR->max[i] >= boundingVolume->max[i]))){
//             return false;
//         }
//     }

//     return true;
// }


//checks to see if the query MBB is outside of the total bounding volume
//Compares the intervals of two hyper-rectangles
bool rangeQueryBoundingVolumeDisjoint(Rectangle * rangeQueryMBR, Rectangle * boundingVolume)
{
    for (unsigned int i=0; i<IDIM; i++)
    {
        if((rangeQueryMBR->max[i] < boundingVolume->min[i]) || (rangeQueryMBR->min[i] > boundingVolume->max[i])){
            return true;
        }
    }

    return false;
}

//Range query search
unsigned int searchRangeQueryWithCandidateSet(Node * root, Rectangle * rangeQueryMBR, 
    Rectangle * boundingVolume, unsigned int depth, Candidates * candSet){

    unsigned int count = 0;

    // Base cases
    //If the root is NULL then the tree is empty or
    //there are no children of the node so return false

    if (root == NULL){
        return 0;
        // printf("\nCase 1: NULL");
    }
    //No overlap with global bounding box
    //i.e., completely disjoint
    //This compares two bounding boxes
    else if(rangeQueryBoundingVolumeDisjoint(rangeQueryMBR, boundingVolume)){
        // printf("\nCase 2: RQ disjoint and outside of bounding volume");
        return 0;
    }
    //Range contains entire cell
    //This compares two bounding boxes
    //This is an optimization for counting range queries 
    //so that traversing to children is unneeded
    //but it does not work when adding points to the candidate set
    /*
    else if(rangeQueryBoundingVolumeContained(rangeQueryMBR, boundingVolume)){
        // printf("\nCase 3: Bounding volume inside RQ");
        updateCandidateSet(root, candSet);
        return root->cnt;
    }
    */
    // range partially overlaps cell
    //This compares a bounding box and a point
    else { 
        // printf("\nCase 4: Partial overlap");
        
        // consider this point
        if (rangeQueryPointContained(rangeQueryMBR, root)){
        // printf("\nCase 5: This point");
        updateCandidateSet(root, candSet);
        // printf("\n[Case 5] Updated candidate set by 1: total: %u", candSet->cnt);
        count++;
        } 
        
        //apply recursively to children
        unsigned int cd = depth % IDIM;

        // in prior code leftPartMBB & rightPartMBB were dynamically allocated
        //to avoid knowing the dimension at compile time
        //But it was much slower due to the recursive mallocs & frees
        Rectangle leftPartMBB = leftpart(cd, root->point, boundingVolume); 
        Rectangle rightPartMBB = rightpart(cd, root->point, boundingVolume); 

        count += searchRangeQueryWithCandidateSet(root->left, rangeQueryMBR, &leftPartMBB, depth+1, candSet);
        count += searchRangeQueryWithCandidateSet(root->right, rangeQueryMBR, &rightPartMBB, depth+1, candSet);
        
        
        
    }

    return count;

}

Rectangle computeBoundingVolume(DTYPE * points, unsigned int numPoints, unsigned int NDIM)
{

    Rectangle globalMBB;

    // printf("\nComputing bounding volume\n");
    //init to first point
    for (unsigned int j=0; j<IDIM; j++)
    {
        globalMBB.min[j]=points[j];
        globalMBB.max[j]=points[j];
    }

    //temp
    // unsigned int idxMin[IDIM]={0,0};
    // unsigned int idxMax[IDIM]={0,0};

    for(unsigned int i=0; i<numPoints; i++)
    {
        for (unsigned int j=0; j<IDIM; j++)
        {
            if(points[(i*NDIM)+j]<globalMBB.min[j]){
                globalMBB.min[j]=points[(i*NDIM)+j];
                // idxMin[j]=i; 
            }
            if(points[(i*NDIM)+j]>globalMBB.max[j]){
                globalMBB.max[j]=points[(i*NDIM)+j];
                // idxMax[j]=i;
            }
        }
    }
    
    printf("\nGlobal bounding volume:");
    for (unsigned int j=0; j<NDIM; j++){
        std::cout<<"\nDim: "<<j<<"min/max: "<<globalMBB.min[j]<<", "<<globalMBB.max[j];
    }

    return globalMBB;

}





//sortDim -- the data dimension to sort tree nodes on
void sortKeyVal(keyValData * keyValDataset, unsigned int NPOINTS, unsigned int sortDim)
{
    //set the sort dimension in each node
    //This will be used in the comparison function
    for (unsigned int i=0; i<NPOINTS; i++){
        keyValDataset[i].sortDim=sortDim;
    }

    //sort sequential
    // qsort(keyValDataset, NPOINTS, sizeof(keyValData), compareKeyValDataStruct);

    //sort parallel
    __gnu_parallel::sort(keyValDataset, keyValDataset+NPOINTS, compareKeyValDataStructParallel);


}




//returns true if the point is contained within the MBR
bool rangeQueryPointContained(Rectangle * rangeQueryMBR, NodeIndices * root)
{
    for (unsigned int i=0; i<IDIM; i++)
    {
        if(!(root->point[i]>=rangeQueryMBR->min[i] && root->point[i]<=rangeQueryMBR->max[i]))
        {
            return false;
        }

    }

    return true;
}


void updateCandidateSet(NodeIndices * root, Candidates * candSet)
{

    if(candSet->cnt==CANDBUFFERSIZE){
        printf("\n[Error] Candidate set overflow! Candidate set size is %u", candSet->cnt);
        exit(0);
    }

    //add buffer overflow logic
    //add current buffer size to the candSet struct as well 
    //(reinventing vectors)
    candSet->candidates[candSet->cnt] = root->pid;
    candSet->cnt++;
}



//Range query search
unsigned int searchRangeQueryWithCandidateSetIdx(int root, NodeIndices * treeNodes, Rectangle * rangeQueryMBR, 
    Rectangle * boundingVolume, unsigned int depth, Candidates * candSet){

    unsigned int count = 0;

    // Base cases
    //If the root is NULL then the tree is empty or
    //there are no children of the node so return false

    if (root == -1){
        return 0;
        // printf("\nCase 1: NULL");
    }
    //No overlap with global bounding box
    //i.e., completely disjoint
    //This compares two bounding boxes
    else if(rangeQueryBoundingVolumeDisjoint(rangeQueryMBR, boundingVolume)){
        // printf("\nCase 2: RQ disjoint and outside of bounding volume");
        return 0;
    }
    //Range contains entire cell
    //This compares two bounding boxes
    //This is an optimization for counting range queries 
    //so that traversing to children is unneeded
    //but it does not work when adding points to the candidate set
    /*
    else if(rangeQueryBoundingVolumeContained(rangeQueryMBR, boundingVolume)){
        // printf("\nCase 3: Bounding volume inside RQ");
        updateCandidateSet(root, candSet);
        return root->cnt;
    }
    */
    // range partially overlaps cell
    //This compares a bounding box and a point
    else { 
        // printf("\nCase 4: Partial overlap");
        
        // consider this point
        if (rangeQueryPointContained(rangeQueryMBR, &treeNodes[root])){
        // printf("\nCase 5: This point");
        updateCandidateSet(&treeNodes[root], candSet);
        // printf("\n[Case 5] Updated candidate set by 1: total: %u", candSet->cnt);
        count++;
        } 
        
        //apply recursively to children
        unsigned int cd = depth % IDIM;

        // in prior code leftPartMBB & rightPartMBB were dynamically allocated
        //to avoid knowing the dimension at compile time
        //But it was much slower due to the recursive mallocs & frees
        Rectangle leftPartMBB = leftpart(cd, treeNodes[root].point, boundingVolume); 
        Rectangle rightPartMBB = rightpart(cd, treeNodes[root].point, boundingVolume); 

        count += searchRangeQueryWithCandidateSetIdx(treeNodes[root].leftIdx, treeNodes, rangeQueryMBR, &leftPartMBB, depth+1, candSet);
        count += searchRangeQueryWithCandidateSetIdx(treeNodes[root].rightIdx, treeNodes, rangeQueryMBR, &rightPartMBB, depth+1, candSet);
        
        
        
    }

    return count;

}

//Similarity search: refine the candidate set using Euclidean distances between
//the candidate set and the query point ID
unsigned int refineCandidateSet(DTYPE * dataset, DTYPE * queryset, unsigned int queryPointId, Candidates * candSet, 
    unsigned int NDIM, DTYPE epsilon, std::vector<unsigned int> *resultSet){
    
    const DTYPE epsilonSq = epsilon*epsilon;

    unsigned int count = 0;
    
    DTYPE tmpdistance = 0;
    unsigned int candidateID;

    //true means that the distance calculation didn't break
    bool breakFlag;
    DTYPE distance = 0;

    for(unsigned int i=0; i<candSet->cnt; i++){
        breakFlag = true;
        distance = 0;
        candidateID = candSet->candidates[i];

        for(unsigned int j=0; j<NDIM; j++){
                tmpdistance = queryset[(queryPointId*NDIM)+j]-dataset[(candidateID*NDIM)+j];
                distance += (tmpdistance*tmpdistance); 
                
                #if SHORTCIRCUIT
                if(distance>epsilonSq){
                    breakFlag = false;
                    break;
                }
                #endif
        }

        if(breakFlag && distance<=epsilonSq){
            count++;

            //If each thread stores key/value pairs in a single vector
            #if COUNTONLY==0 && RESULTSETSTORAGECPU==1
            resultSet->push_back(queryPointId);
            resultSet->push_back(candidateID);
            #endif

            //If we store the results in the vector where we have one vector per query point
            #if COUNTONLY==0 && RESULTSETSTORAGECPU==2
            resultSet->push_back(candidateID);
            #endif

        }
    }

    return count;
}

void initPointsInPreallocatedNodes(NodeIndices * treeNodes, DTYPE * dataset, unsigned int NPOINTS, unsigned int NDIM)
{
    for(unsigned int i=0; i<NPOINTS; i++){
        treeNodes[i].pid = i;
        treeNodes[i].point = &dataset[i*NDIM];
    }
}


//This does not start at the root each insertation so it cannot be used for 
//counting range queries
//Note: Because this does not start at the root node, it must use the exact median each time it splits
//or some points will not get indexed

//treeNodes is the base of the array, so it doesn't change with the recursive call
//return type (int) is the index into treeNodes
int insertAllWithIdx(int root, NodeIndices * treeNodes, keyValData * datasetKeyVal, unsigned int NPOINTS, unsigned int depth, unsigned int * maxDepth, unsigned int * insertedNodeDepth, unsigned int NDIM, bool allowDuplicates)
{
    
    // printf("\nNPOINTS in function: %u", NPOINTS);

    if(NPOINTS==0){
        return root;
    }

    //used to determine the depth of the node
    *insertedNodeDepth = depth;

    //keep track of max depth for statistics on tree construction
    if(*maxDepth<*insertedNodeDepth){
        *maxDepth = *insertedNodeDepth;
    }

    //compute the dimension:
    unsigned int currentDim = depth % IDIM;

    //Sort by dimension
    sortKeyVal(datasetKeyVal, NPOINTS, currentDim);

    unsigned int mid = NPOINTS/2;
    DTYPE * point = datasetKeyVal[mid].point;

    //initialize node and set the root
    initNode(&treeNodes[datasetKeyVal[mid].pid]);
    root = datasetKeyVal[mid].pid;

    //If allowDuplicates is false and
    //if there are duplicate points then we throw an error and exit
    if(allowDuplicates==false){
        // if(checkEqualPoints(point, root->point)){
        if(checkEqualPoints(point, datasetKeyVal[mid].point)){
            printf("\nThe two points are equal and the kd-tree cannot use duplicates\n");
            printf("\nPoint id being inserted: %u", datasetKeyVal[mid].pid);
            printf("\nPoint id already in the tree: %u", treeNodes[datasetKeyVal[mid].pid].pid);
            exit(0);
        }
    }
    
    unsigned int NPOINTSLEFT = mid;
    unsigned int NPOINTSRIGHT = NPOINTS - (mid+1);
    
    
    treeNodes[root].leftIdx = insertAllWithIdx(treeNodes[root].leftIdx, treeNodes, datasetKeyVal, NPOINTSLEFT, depth+1, maxDepth, insertedNodeDepth, NDIM, allowDuplicates);
    treeNodes[root].rightIdx = insertAllWithIdx(treeNodes[root].rightIdx, treeNodes, &datasetKeyVal[mid+1], NPOINTSRIGHT, depth+1, maxDepth, insertedNodeDepth, NDIM, allowDuplicates);  
    //update the counter so that this node has a count of all subtrees
    //This will not provide an accurate count because we do not insert at the node each time
    treeNodes[root].cnt++;
    return root;
}



uint64_t searchAndRefineWithMode(const unsigned int MODE, Node * root, int rootIdx, NodeIndices * treeNodes, DTYPE * dataset, DTYPE * queryset, 
    unsigned int NDIM, unsigned int NPOINTS, unsigned int NQPOINTS, DTYPE epsilon, void * resultSetVoid
    // std::vector<unsigned int> *resultSet
    )
{

    //cast the void type to either a vector or unsigned int
    #if COUNTONLY==0
    std::vector<unsigned int> *resultSet = (std::vector<unsigned int> *)resultSetVoid;
    #endif

    #if COUNTONLY==1
    unsigned int * resultSet = (unsigned int*)resultSetVoid;
    #endif

    //populate globalMBB which is the global bounding volume
    // Rectangle globalMBB = computeBoundingVolume(queryset, NQPOINTS, NDIM);
    Rectangle globalMBB = computeBoundingVolume(dataset, NPOINTS, NDIM);


    //Candidate set buffer to be refined
    //In the parallel version we will need one buffer per thread (NUMCPUTHREADS total threads)
    Candidates * candSet = (Candidates*)malloc(sizeof(Candidates)*NUMCPUTHREADS);
    Candidates * candSetPtr;
    

    //Reuse the queryMBB for the searches across all points
    //Need one queryMBB for each point when running in parallel
    Rectangle queryMBB[NUMCPUTHREADS];
    Rectangle * queryMBBPtr;

    //Various counters
    unsigned int maxCnt = 0;
    uint64_t totalCountSimSearch = 0;
    uint64_t totalCountRQ = 0;
    uint64_t totalCountRQCandSet = 0; //Sanity check that we haven't missed any points
    unsigned int startingDepth=0;
    unsigned int countRQ=0;

    double totalTimeRefineCandidates = 0;


    


    
    #pragma omp parallel for num_threads(NUMCPUTHREADS) private(candSetPtr, queryMBBPtr, countRQ, startingDepth)\
    reduction(+:totalCountSimSearch, totalCountRQ, totalCountRQCandSet, totalTimeRefineCandidates) reduction(max: maxCnt)\
    schedule(guided)
    for (unsigned int i=0; i<NQPOINTS; i++){

        unsigned int tid = omp_get_thread_num();




        candSetPtr = &candSet[tid];
        queryMBBPtr = &queryMBB[tid];

        //reset candidate set counter
        candSetPtr->cnt = 0;
        

        countRQ=0;
        generateQueryMBBFromPoint(queryMBBPtr, queryset+(i*NDIM), epsilon);   
        
        //When using Node struct
        if(MODE==0){
            //count from the range query 
            countRQ = searchRangeQueryWithCandidateSet(root, queryMBBPtr, &globalMBB, 
            startingDepth, candSetPtr);
        }
        //When using NodeIndices struct
        if(MODE==1){
            //count from the range query 
            countRQ = searchRangeQueryWithCandidateSetIdx(rootIdx, treeNodes, queryMBBPtr, &globalMBB, 
            startingDepth, candSetPtr);
        }

        //This should be the same for the search mode that can do counting range queries
        totalCountRQ+= (uint64_t)countRQ;
        totalCountRQCandSet+= (uint64_t)candSetPtr->cnt;

        //count after refining for similarity search (using Euclidean distance instead of range query)
        double tstartRefine=omp_get_wtime();

        //if we do not store the result set
        #if COUNTONLY==0 && RESULTSETSTORAGECPU==0 
        unsigned int countSimSearch = refineCandidateSet(dataset, queryset, i, candSetPtr, NDIM, epsilon, NULL);
        #endif

        //if we use one vector per thread
        #if COUNTONLY==0 && RESULTSETSTORAGECPU==1
        unsigned int countSimSearch = refineCandidateSet(dataset, queryset, i, candSetPtr, NDIM, epsilon, &resultSet[tid]);
        #endif

        //if we use one vector per query point
        #if COUNTONLY==0 && RESULTSETSTORAGECPU==2
        unsigned int countSimSearch = refineCandidateSet(dataset, queryset, i, candSetPtr, NDIM, epsilon, &resultSet[i]);
        #endif

        //if we do counting range queries
        #if COUNTONLY==1 
        unsigned int countSimSearch = refineCandidateSet(dataset, queryset, i, candSetPtr, NDIM, epsilon, NULL);
        resultSet[i] = countSimSearch;
        #endif
        
        totalCountSimSearch+= (uint64_t)countSimSearch;
        
        double tendRefine=omp_get_wtime();
        totalTimeRefineCandidates+= tendRefine - tstartRefine;

        
        //For testing epsilon and size of the candidate set
        if(maxCnt<countRQ){
            maxCnt = countRQ;
        }       


        if(i%10000==0){
        printf("\nPercentage complete: %f", ((i*1.0/NQPOINTS)*100.0));
        }

    }

    printf("\nTotal count range query: %lu", totalCountRQ);
    printf("\nTotal count range query (candidate set; should be same as above): %lu", totalCountRQCandSet);
    printf("\nTotal count similarity search: %lu", totalCountSimSearch);
    printf("\nMaximum count in the candidate set: %u", maxCnt);

    //this only makes sense if this is sequential
    if(NUMCPUTHREADS==1){
        printf("\nTotal Time refine candidate set (only, excluding searching the index): %f", totalTimeRefineCandidates);
    }


    return totalCountSimSearch;
    

}

