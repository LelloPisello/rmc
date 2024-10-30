#ifndef RMC_h_
#define RMC_h_
#include <stdint.h>


//i'm sorry this needed to be quick i didnt have time to make a proper debug system
#ifdef DEBUG
    #include <stdio.h>
    #define DEBUGPRINT(...) fprintf(stderr, __VA_ARGS__)
#else
    #define DEBUGPRINT(...) 
#endif

#ifdef RMC_USE_FLOAT 
#define RmcFloat float
#else
#define RmcFloat double
#endif

#define RMC_DEFINE_HANDLE(handlename) typedef struct handlename##_s *handlename;

//Result type, used for every function that can return an error
typedef enum {
    RMC_SUCCESS = 0,
    
    //TODO: make errors flag their corresponding bit
    RMC_ERROR_GENERIC_BIT = 1 << 31,
    RMC_ERROR_TENSOR_BIT = 1 << 30,

    //generic failures
    RMC_ERROR_GENERIC_EMPTY_HANDLE,
    RMC_ERROR_GENERIC_NULLPTR,
    RMC_ERROR_GENERIC_INVALID_PTR,
    RMC_ERROR_GENERIC_FAILED_ALLOCATION,
    RMC_ERROR_GENERIC_ZERO_SIZE,

    //manifold specific errors
    RMC_ERROR_MANIFOLD_METRIC_FIELD_CALCULATION_FAILURE,
    RMC_ERROR_MANIFOLD_CHRISTOFFEL_FIELD_CALCULATION_FAILURE,

    //tensor operation specific errors

    //covector + vector
    RMC_ERROR_TENSOR_TYPES_INCOMPATIBLE,

    //attempting to lower index on a covector or viceversa
    RMC_ERROR_TENSOR_TYPE_WRONG,
    

} RmcError;

RMC_DEFINE_HANDLE(RmcManifold)

#undef RMC_DEFINE_HANDLE

#endif