#ifndef RMC_MANIFOLD_h_
#define RMC_MANIFOLD_h_

#include "geo.h"
#include "rmc.h"

//MANIFOLD:
/*
Object containing data about the manifold (scene)
*/

typedef RmcFloat RmcMetricOutput[3][3];

/*
    pointer to a generating function for a variable metric.
    will be evaluated for each cell in the 3D buffer and saved.
*/
typedef void (*RmcMetricGenerator)(const RmcCoordinates* coords, RmcMetricOutput* output);

typedef struct {
    //coordinate bounds [-fBounds, fBounds]
    RmcFloat fBounds;

    /*
    every field on the manifold (metric, christoffel symbols)
    is saved onto a 3D {uRes, uRes, uRes} buffer
    */
    uint32_t uResolution;

    //pointer to the chosen generating function for the metric
    RmcMetricGenerator fpMetricGenerator;
} RmcManifoldCreateInfo;

//creates a manifold based on the data contained in the create info structure
RmcError rmcManifoldCreate(const RmcManifoldCreateInfo* manifoldInfo, RmcManifold* result);

//obtains length of vector
RmcError rmcManifoldTensorGetLength(RmcManifold manifold, const RmcCoordinates* coords, const RmcTensor* vector, RmcFloat* result);

//vector->covector
RmcError rmcManifoldTensorLowerIndex(RmcManifold manifold, const RmcCoordinates* coords, const RmcTensor* vector, RmcTensor* result);
//covector->vector
RmcError rmcManifoldTensorRaiseIndex(RmcManifold manifold, const RmcCoordinates* coords, const RmcTensor* vector, RmcTensor* result);

//frees the manifold and its associated buffers
void rmcManifoldDestroy(RmcManifold manifold);

#endif