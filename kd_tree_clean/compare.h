// #include <cstdint>
// #include "params.h"
// #include "structs.h"


//Comparison functions
int compare(const void* a, const void* b) {
   return (*(int*)a - *(int*)b);
}

//sort descending
int compareByDimVariance(const void* a, const void* b)
{
    double l = ((struct dim_reorder_sort*)a)->variance;
    double r = ((struct dim_reorder_sort*)b)->variance;

    if (l > r){
      return -1;
   }
   else if (l < r){
      return 1;
   }
   else{
      return 0;
   }
 
}

//sort key/value struct on the sorting dimension
int compareKeyValDataStruct(const void* a, const void* b) {
   
   unsigned int sortDimA = ((struct keyValData*)a)->sortDim; 
   unsigned int sortDimB = ((struct keyValData*)b)->sortDim; 
   double l = ((struct keyValData*)a)->point[sortDimA]; 
   double r = ((struct keyValData*)b)->point[sortDimB]; 
   
   if (l < r){
      return -1;
   }
   else if (l > r){
      return 1;
   }
   else{
      return 0;
   }
}



bool compareDescendingKeyVal(struct keyValPair l, struct keyValPair r)
{
   return l.key > r.key;
}


//sort key/value struct on the sorting dimension for parallel sort
bool compareKeyValDataStructParallel(struct keyValData a, struct keyValData b) {
   
   unsigned int sortDimA = a.sortDim; 
   unsigned int sortDimB = b.sortDim; 
   double l = a.point[sortDimA]; 
   double r = b.point[sortDimB];  
   
   return (l < r);
}

//sort key/value struct on the first dimensional coordinate
bool compareKeyValDataStructCompareOnOneCoord(struct keyVal a, struct keyVal b) {


   
   DTYPE l = a.pointSortedDim; 
   DTYPE r = b.pointSortedDim;
   
   return l < r;   

}






//sort key/value struct on the first dimensional coordinate
bool compareKeyValDataStructLowerUpperBound(struct keyVal &a, struct keyVal &b) {


   
   DTYPE l = a.pointSortedDim; 
   DTYPE r = b.pointSortedDim;
   
   return l < r;   

}


//sort key/value result set pairs
bool compareKeyValResult(struct keyValResultSet &a, struct keyValResultSet &b){
   return a.key < b.key;   
}
 


