#ifndef RMC_GEO_h_
#define RMC_GEO_h_

#include "rmc.h"


//TENSORS

//defining length and products without a metric is impossible
//reserved for later implementation

//Simple enum for vector-covector distinction
typedef enum {
    RMC_TENSOR_TYPE_VECTOR,
    RMC_TENSOR_TYPE_COVECTOR
} RmcTensorType;

//a 3-tuple meant to represent vector components/coordinates
typedef union {
    RmcFloat u[3];
    struct { RmcFloat x, y, z; };
} RmcCoordinates;

//3d (t, !t)tensor (covector or vector)
typedef struct {
    RmcTensorType type;
    RmcCoordinates components;
} RmcTensor;

//sum between two tensor of the same type
RmcError rmcTensorSum(const RmcTensor* a, const RmcTensor* b, RmcTensor* result);
//mul between scalar and tensor
RmcError rmcTensorMul(RmcFloat r, const RmcTensor* t, RmcTensor* result);

#endif