#include "params.h"
#include <stdint.h>
#include <limits.h>
#include <stdio.h>
// #include <dpu>


#ifndef __COMMON_H__
#define __COMMON_H__

// #define XSTR(x) STR(x)
// #define STR(x) #x

//DPU variables for host
#define DPU_BUFFER dpu_mram_buffer
#define DPU_CACHES dpu_wram_caches
#define DPU_RESULTS dpu_wram_results



#define BUFFER_SIZE 16

//Structure used by both the host and the dpu to communicate information
//the result (a sum and the number of cycles used by each tasklet)
typedef struct {
    uint64_t sum;
    uint64_t cycles;
} dpu_result_t;

typedef struct {
    dpu_result_t tasklet_result[NR_TASKLETS];
} dpu_results_t;

//In arguments used across a DPU set
typedef struct{
    uint64_t dpu_id;
    uint64_t global_dpu_id; //for debugging
    uint64_t num_dpus;
    uint64_t num_dim;
    uint64_t sort_dim;
    uint64_t numPointsDatasetPartition;
    uint64_t numPointsQuerysetPartition;
    DTYPE epsilon;  
} dpu_arguments_t;





#endif
