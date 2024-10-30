#include "manifold.h"
#include "geo.h"
#include "rmc.h"
#include <stdlib.h>
#include <math.h>


/*
    QUICK COMMENT ABOUT MEMORY USAGE

    the choice was either huge computation time or huge memory print.
    i went with the second.
    I could either compute the christoffel symbols every time they were needed or 
    I could store them in a huge buffer.
    The first would have made the program run incredibly slowly.
*/

//self explanatory. first index is upper index, second and third lower left to right
/*
    ._ a
    |  bc

    [a][b][c]
*/
typedef RmcFloat RmcChristoffelSymbol[3][3][3];

struct RmcManifold_s {
    RmcFloat fBounds;
    uint32_t uResolution;
    uint32_t uSize;//entries in buffers

    //for lowering
    RmcMetricOutput *metricField;
    //for raising
    RmcMetricOutput *inverseMetricField;
    RmcChristoffelSymbol *christoffelField;

    

};

static void _getMetricAt(RmcManifold manifold, const RmcCoordinates* coords, RmcMetricOutput* output) {
    const float dXi2 = manifold->fBounds / manifold->uResolution;
    RmcCoordinates adjustedCoords = {
        .x = (coords->x + manifold->fBounds - dXi2) / manifold->fBounds / 2 * manifold->uResolution,
        .y = (coords->y + manifold->fBounds - dXi2) / manifold->fBounds / 2 * manifold->uResolution,
        .z = (coords->z + manifold->fBounds - dXi2) / manifold->fBounds / 2 * manifold->uResolution,
    };  
    uint32_t i = adjustedCoords.x;
    i += adjustedCoords.y * manifold->uResolution;
    i += adjustedCoords.z * manifold->uResolution * manifold->uResolution;

    for(uint32_t j = 0; j < 3; ++j) {
        for(uint32_t k = 0; k < 3; ++k) {
            (*output)[j][k] = manifold->metricField[i][j][k];
        }
    }
}

//calculate values for metric in buffer
static RmcError _fillMetric(RmcManifold manifold, RmcMetricGenerator func) {
    //width of cube 
    const RmcFloat dXi = (manifold->fBounds * 2) / manifold->uResolution;

    for(uint32_t i = 0; i < manifold->uSize; ++i) {
        RmcCoordinates coords = {
            .x = (i % manifold->uResolution),
            .y = (i % (manifold->uResolution * manifold->uResolution)) / manifold->uResolution,
            .z = (i / (manifold->uResolution * manifold->uResolution))
        };
        //get coords
        coords.x *= dXi;
        coords.y *= dXi;
        coords.z *= dXi;
        //center coords
        const RmcFloat offset = dXi * 0.5 * (manifold->uResolution + 1);
        coords.x -= offset;
        coords.y -= offset;
        coords.z -= offset;

        //fill metric in spot
        func(&coords, &manifold->metricField[i]);
    }
    return RMC_SUCCESS;
}

RmcError rmcManifoldCreate(const RmcManifoldCreateInfo* manifoldInfo, RmcManifold* result){
    if(!manifoldInfo || !manifoldInfo->fpMetricGenerator) {
        return RMC_ERROR_GENERIC_NULLPTR;
    }

    if(!manifoldInfo->uResolution) {
        return RMC_ERROR_GENERIC_ZERO_SIZE;
    }
#define ALIAS (*result)

    DEBUGPRINT("Creating manifold...\n");

    *result = malloc(sizeof(**result));

    ALIAS->fBounds = manifoldInfo->fBounds;
    ALIAS->uResolution = manifoldInfo->uResolution;

    ALIAS->uSize = 
        ALIAS->uResolution *
        ALIAS->uResolution *
        ALIAS->uResolution;

    //create buffers

    DEBUGPRINT("\tAllocating metric field for manifold...\n");
    ALIAS->metricField = malloc(
        ALIAS->uSize *
        sizeof(RmcMetricOutput)
    );

    if(!ALIAS->metricField) {
        DEBUGPRINT("\tFailed to allocate metric field\n\n");
        free(ALIAS);
        return RMC_ERROR_GENERIC_FAILED_ALLOCATION;
    }
    DEBUGPRINT("\tOK\n\n");

    DEBUGPRINT("\tAllocating christoffel field for manifold...\n");
    ALIAS->christoffelField = malloc(
        ALIAS->uSize *
        sizeof(RmcChristoffelSymbol)
    );

    if(!ALIAS->christoffelField) {
        DEBUGPRINT("\tFailed to allocate christoffel field\n\n\t");
        free(ALIAS->metricField);
        free(ALIAS);
        return RMC_ERROR_GENERIC_FAILED_ALLOCATION;
    }
    DEBUGPRINT("\tOK\n\n\t");

    DEBUGPRINT("Filling metric using the generator...\n\t");
    RmcError temp = _fillMetric(ALIAS, manifoldInfo->fpMetricGenerator);

    if(temp) {
        DEBUGPRINT("Failed to fill manifold using generator");
        return temp;
    }
    DEBUGPRINT("OK\n\nSuccessfully created manifold\n\n");

#undef ALIAS

    return RMC_SUCCESS;
}

void rmcManifoldDestroy(RmcManifold manifold) {
    free(manifold->christoffelField);
    free(manifold->metricField);
    free(manifold);
    DEBUGPRINT("Destroyed a manifold\n");
}

RmcError rmcManifoldTensorLowerIndex(RmcManifold manifold, const RmcCoordinates* coords, const RmcTensor* vector, RmcTensor* result) {

    return RMC_SUCCESS;
}

RmcError rmcManifoldTensorGetLength(RmcManifold manifold, const RmcCoordinates* coords, const RmcTensor* vector, RmcFloat* result) {
    if(!manifold) {
        return RMC_ERROR_GENERIC_EMPTY_HANDLE;
    }
    if(!coords || !vector || !result) {
        return RMC_ERROR_GENERIC_NULLPTR;
    }
    if(vector->type == RMC_TENSOR_TYPE_COVECTOR) {
        return RMC_ERROR_TENSOR_TYPE_WRONG;
    }

    RmcMetricOutput metric;
    _getMetricAt(manifold, coords, &metric);
    RmcFloat length = 0;
    for(uint32_t i = 0; i < 3; ++i) {
        for(uint32_t j = 0; j < 3; ++j) {
            length += metric[i][j] * vector->components.u[i] * vector->components.u[j];
        }
    }
    *result = sqrt(length);
    return RMC_SUCCESS;
}