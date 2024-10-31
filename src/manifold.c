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
    uint64_t uMetricResolution;
    uint64_t uMetricSize;//entries in buffers

    //metricresolution + 1
    uint64_t uChristoffelResolution;
    //size of christoffel buffer (+1 for each dimension)
    uint64_t uChristoffelSize;


    //for lowering
    RmcMetricOutput *metricField;
    //for raising
    RmcMetricOutput *inverseMetricField;
    RmcChristoffelSymbol *christoffelField;

    

};

//i have copied this without remorse from stackoverflow
static void _invertMetric(const RmcMetricOutput* metric, RmcMetricOutput* result) {
#define A (*metric)
#define R (*result)
    const RmcFloat det = +A[0][0]*(A[1][1]*A[2][2]-A[2][1]*A[1][2])
                        -A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
                        +A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]);
    const RmcFloat invdet = 1. / det;


    R[0][0] =  (A[1][1]*A[2][2]-A[2][1]*A[1][2])*invdet;
    R[1][0] = -(A[0][1]*A[2][2]-A[0][2]*A[2][1])*invdet;
    R[2][0] =  (A[0][1]*A[1][2]-A[0][2]*A[1][1])*invdet;
    R[0][1] = -(A[1][0]*A[2][2]-A[1][2]*A[2][0])*invdet;
    R[1][1] =  (A[0][0]*A[2][2]-A[0][2]*A[2][0])*invdet;
    R[2][1] = -(A[0][0]*A[1][2]-A[1][0]*A[0][2])*invdet;
    R[0][2] =  (A[1][0]*A[2][1]-A[2][0]*A[1][1])*invdet;
    R[1][2] = -(A[0][0]*A[2][1]-A[2][0]*A[0][1])*invdet;
    R[2][2] =  (A[0][0]*A[1][1]-A[1][0]*A[0][1])*invdet;
#undef R
#undef A
}

static void _getChristoffelAt(RmcManifold manifold, const RmcCoordinates* coords, RmcChristoffelSymbol* result) {
    for(uint32_t i = 0; i < 27; ++i) {
        (*result)[i % 3][i / 3 % 3][i / 9] = 0.0;
    }
    
    const RmcFloat r = (manifold->uMetricResolution / 2.0);
    const RmcCoordinates fIndex = {
        .x = r * (coords->x / manifold->fBounds / (1 + 1.0 / manifold->uMetricResolution) + 1),
        .y = r * (coords->y / manifold->fBounds / (1 + 1.0 / manifold->uMetricResolution) + 1),    
        .z = r * (coords->z / manifold->fBounds / (1 + 1.0 / manifold->uMetricResolution) + 1)
    };
    RmcCoordinates delta = {
        .x = fIndex.x - (uint32_t)fIndex.x,
        .y = fIndex.y - (uint32_t)fIndex.y,
        .z = fIndex.z - (uint32_t)fIndex.z
    };

    const uint64_t 
        ix = fIndex.x,
        iy = fIndex.y,
        iz = fIndex.z;

#define CHRISTOFFELINDEX(x, y, z)(x + y * manifold->uChristoffelResolution + z * manifold->uChristoffelResolution * manifold->uChristoffelResolution)

    for(uint32_t d = 0; d < 8; ++d) {
        const uint32_t 
            onx = d & 1 ? 1 : 0,
            ony = d & 2 ? 1 : 0,
            onz = d & 4 ? 1 : 0;
        for(uint32_t i = 0; i < 27; ++i) {
            (*result)[i % 3][i / 3 % 3][i / 9] += 
                (onx ? delta.x : 1 - delta.x) *
                (ony ? delta.y : 1 - delta.y) *
                (onz ? delta.z : 1 - delta.z) *
                manifold->christoffelField[CHRISTOFFELINDEX(
                    ix + onx, 
                    iy + ony, 
                    iz + onz
                )][i % 3][i / 3 % 3][i / 9];

        }
        
    }

#undef CHRISTOFFELINDEX
}

static void _getMetricAt(RmcManifold manifold, RmcMetricOutput* metricField, const RmcCoordinates* coords, RmcMetricOutput* output) {
    const RmcFloat r = (manifold->uMetricResolution - 1) / 2.0 / manifold->fBounds;

    for(uint32_t i = 0; i < 3; ++i) {
        for(uint32_t j = 0; j < 3; ++j) {
            (*output)[i][j] = 0.0;
        }
    }    

    const RmcCoordinates fIndex = { 
        .x = r * (coords->x + manifold->fBounds),
        .y = r * (coords->y + manifold->fBounds),
        .z = r * (coords->z + manifold->fBounds)
    };
    const RmcCoordinates delta = {
        .x = fIndex.x - (uint32_t)fIndex.x,
        .y = fIndex.y - (uint32_t)fIndex.y,
        .z = fIndex.z - (uint32_t)fIndex.z
    };

    const uint64_t 
        ix = fIndex.x,
        iy = fIndex.y,
        iz = fIndex.z;

#define METRICINDEX(x, y, z) (x + y * manifold->uMetricResolution + z * manifold->uMetricResolution * manifold->uMetricResolution)

    //kill me I am not writing everything down or using binary digits
    for(uint32_t d = 0; d < 8; ++d) {
        //d & 1, d & 2, d & 4
        //x      y      z
        const uint32_t 
            onx = d & 1 ? 1 : 0,
            ony = d & 2 ? 1 : 0,
            onz = d & 4 ? 1 : 0;
        for(uint32_t i = 0; i < 3; ++i) {
            for(uint32_t j = 0; j < 3; ++j) {
                (*output)[i][j] += 
                (onx ? delta.x : 1 - delta.x) *
                (ony ? delta.y : 1 - delta.y) *
                (onz ? delta.z : 1 - delta.z) *
                metricField
                [METRICINDEX(
                    ix + onx, 
                    iy + ony, 
                    iz + onz)][i][j];
            }
        }
    }

#undef METRICINDEX
}


//calculate values for metric in buffer
static void _fillMetric(RmcManifold manifold, RmcMetricGenerator func) {
    //width of cube 
    const RmcFloat r = 2 * (manifold->uMetricResolution - 1) * manifold->fBounds;

    for(uint64_t i = 0; i < manifold->uMetricSize; ++i) {
        RmcCoordinates coords = {
            .x = (i % manifold->uMetricResolution),
            .y = (i / manifold->uMetricResolution % manifold->uMetricResolution),
            .z = (i / manifold->uMetricResolution / manifold->uMetricResolution)
        };
        //get coords
        coords.x *= r;
        coords.y *= r;
        coords.z *= r;
        //center coords
        coords.x -= manifold->fBounds;
        coords.y -= manifold->fBounds;
        coords.z -= manifold->fBounds;

        //fill metric in spot
        func(&coords, &manifold->metricField[i]);
        _invertMetric(&manifold->metricField[i], &manifold->inverseMetricField[i]);
    }
}

//for levi-civita
static void _fillChristoffel(RmcManifold manifold) {
    const RmcFloat r = manifold->fBounds * (1.0 + 1.0 / manifold->uMetricResolution);

    for(uint64_t i = 0; i < manifold->uChristoffelSize; ++i) {
        const uint64_t 
            ix = i % manifold->uChristoffelResolution,
            iy = i / manifold->uChristoffelResolution % manifold->uChristoffelResolution,
            iz = i / manifold->uChristoffelResolution / manifold->uChristoffelResolution;
        if(
            !ix || !iy || !iz ||
            ix == manifold->uMetricResolution ||
            iy == manifold->uMetricResolution || 
            iz == manifold->uMetricResolution
        ) {
            //I love nested fors i love nested fors i love nested fors
            //the compiler will most definitely optimize this out 
            for(uint32_t j = 0; j < 27; ++j) {
                manifold->christoffelField[i][j % 3][j / 3 % 3][j / 9] = 0.0;
            }
        } else {
            RmcCoordinates baseCoords;
            
            //half of width of brick
            const RmcFloat dOffset = manifold->fBounds / manifold->uMetricResolution;

            {
                const RmcFloat r = (manifold->fBounds + manifold->fBounds / manifold->uMetricResolution);
                baseCoords = (RmcCoordinates){
                    .x = r * (2.0 * ix / manifold->uMetricResolution - 1),
                    .y = r * (2.0 * iy / manifold->uMetricResolution - 1),
                    .z = r * (2.0 * iz / manifold->uMetricResolution - 1)
                };
            }

            //metric tensor derivatives
            RmcMetricOutput metricDerivatives[3];
            {
                RmcMetricOutput temp;
                for(uint32_t j = 0; j < 3; ++j) {
                    baseCoords.u[j] += dOffset;
                    _getMetricAt(manifold, manifold->metricField, &baseCoords, &temp);
                    for(uint32_t k = 0; k < 9; ++k) {
                        metricDerivatives[j][k % 3][k / 3] += temp[k % 3][k / 3];
                    }

                    baseCoords.u[j] -= 2 * dOffset;
                    _getMetricAt(manifold, manifold->metricField, &baseCoords, &temp);
                    baseCoords.u[j] += dOffset;
                    for(uint32_t k = 0; k < 9; ++k) {
                        metricDerivatives[j][k % 3][k / 3] -= temp[k % 3][k / 3];
                        metricDerivatives[j][k % 3][k / 3] *= 1.0 / (dOffset * 2);
                    }
                }
            }

            RmcMetricOutput inverseMetric;
            _getMetricAt(manifold, manifold->inverseMetricField, &baseCoords, &inverseMetric);
            
            //i couldnt resist i hate tensor calculus
            for(uint32_t j = 0; j < 27; ++j) {
                manifold->christoffelField[i][j % 3][j / 3 % 3][j / 9] = 0.0;
                for(uint32_t m = 0; m < 3; ++m) {
                    
                    manifold->christoffelField[i][j % 3][j / 3 % 3][j / 9] +=
                    0.5 * inverseMetric[j % 3][m] * (
                        metricDerivatives[j / 3 % 3][m][j / 9] +
                        metricDerivatives[j / 9][m][j / 3 % 3] -
                        metricDerivatives[m][j / 3 % 3][j / 9]
                    );

                    //manifold->christoffelField[i][j][k][l] = 
                }
            }
        }
    }
}

RmcError rmcManifoldCreate(const RmcManifoldCreateInfo* manifoldInfo, RmcManifold* result){
    if(!manifoldInfo || !manifoldInfo->fpMetricGenerator) {
        return RMC_ERROR_GENERIC_NULLPTR;
    }

    if(!manifoldInfo->uMetricResolution || manifoldInfo->fBounds == 0.0) {
        return RMC_ERROR_GENERIC_ZERO_SIZE;
    }

    if(manifoldInfo->fBounds < 0) {
        return RMC_ERROR_GENERIC_NEGATIVE_SIZE;
    }
#define ALIAS (*result)

    DEBUGPRINT("Creating manifold...\n");

    *result = malloc(sizeof(**result));

    ALIAS->fBounds = manifoldInfo->fBounds;
    ALIAS->uMetricResolution = manifoldInfo->uMetricResolution;
    ALIAS->uChristoffelResolution = ALIAS->uMetricResolution + 1;

    ALIAS->uMetricSize = 
        ALIAS->uMetricResolution *
        ALIAS->uMetricResolution *
        ALIAS->uMetricResolution;

    //(uMetricResolution + 1) ^ 3
    //i know this is horrible but i don't care
    ALIAS->uChristoffelSize = 
        ALIAS->uChristoffelResolution *
        ALIAS->uChristoffelResolution *
        ALIAS->uChristoffelResolution;

    DEBUGPRINT("\tLength of metric buffer: %lu\n",
    (long unsigned)ALIAS->uMetricSize);

    //create buffers

    DEBUGPRINT("\tAllocating metric field for manifold...\n");
    ALIAS->metricField = malloc(
        ALIAS->uMetricSize *
        sizeof(RmcMetricOutput)
    );

    if(!ALIAS->metricField) {
        DEBUGPRINT("\tFailed to allocate metric field\n\n");
        free(ALIAS);
        return RMC_ERROR_GENERIC_FAILED_ALLOCATION;
    }
    DEBUGPRINT("\tOK\n\n");

    DEBUGPRINT("\tAllocating inverse metric field for manifold...\n");
    ALIAS->inverseMetricField = malloc(
        ALIAS->uMetricSize *
        sizeof(RmcMetricOutput)
    );

    if(!ALIAS->inverseMetricField) {
        DEBUGPRINT("\tFailed to allocate inverse metric field\n\n");
        free(ALIAS->metricField);
        free(ALIAS);
        return RMC_ERROR_GENERIC_FAILED_ALLOCATION;
    }
    DEBUGPRINT("\tOK\n\n");

    DEBUGPRINT("\tAllocating christoffel field for manifold...\n");
    ALIAS->christoffelField = malloc(
        ALIAS->uChristoffelSize *
        sizeof(RmcChristoffelSymbol)
    );

    if(!ALIAS->christoffelField) {
        DEBUGPRINT("\tFailed to allocate christoffel field\n\n\t");
        free(ALIAS->metricField);
        free(ALIAS->inverseMetricField);
        free(ALIAS);
        return RMC_ERROR_GENERIC_FAILED_ALLOCATION;
    }
    DEBUGPRINT("\tOK\n\n\t");

    DEBUGPRINT("Filling metric using the generator...\n\t");
    _fillMetric(ALIAS, manifoldInfo->fpMetricGenerator);

    DEBUGPRINT("OK\n\n\tFilling christoffel field using the metric\n\t");
    _fillChristoffel(ALIAS);
    DEBUGPRINT("OK\n\nSuccessfully created manifold\n\n");

#undef ALIAS

    return RMC_SUCCESS;
}

void rmcManifoldDestroy(RmcManifold manifold) {
    free(manifold->christoffelField);
    free(manifold->metricField);
    free(manifold->inverseMetricField);
    free(manifold);
    DEBUGPRINT("Destroyed a manifold\n");
}

RmcError rmcManifoldTensorLowerIndex(RmcManifold manifold, const RmcCoordinates* coords, const RmcTensor* vector, RmcTensor* result) {
    if(!manifold) {
        return RMC_ERROR_GENERIC_EMPTY_HANDLE;
    }
    if(!coords || !vector || !result) {
        return RMC_ERROR_GENERIC_NULLPTR;
    }
    if(vector->type != RMC_TENSOR_TYPE_VECTOR) {
        return RMC_ERROR_TENSOR_TYPE_WRONG;
    }

    RmcMetricOutput metric;
    _getMetricAt(manifold, manifold->metricField, coords, &metric);
    result->type = RMC_TENSOR_TYPE_COVECTOR;

    for(uint64_t i = 0; i < 3; ++i) {
        result->components.u[i] = 0;
        for(uint64_t j = 0; j < 3; ++j) {
            result->components.u[i] += metric[i][j] * vector->components.u[j];
        }
    }

    return RMC_SUCCESS;
}

RmcError rmcManifoldTensorRaiseIndex(RmcManifold manifold, const RmcCoordinates* coords, const RmcTensor* vector, RmcTensor* result) {
    if(!manifold) {
        return RMC_ERROR_GENERIC_EMPTY_HANDLE;
    }
    if(!coords || !vector || !result) {
        return RMC_ERROR_GENERIC_NULLPTR;
    }
    if(vector->type != RMC_TENSOR_TYPE_COVECTOR) {
        return RMC_ERROR_TENSOR_TYPE_WRONG;
    }

    RmcMetricOutput metric;
    _getMetricAt(manifold, manifold->inverseMetricField, coords, &metric);
    result->type = RMC_TENSOR_TYPE_VECTOR;

    for(uint64_t i = 0; i < 3; ++i) {
        result->components.u[i] = 0;
        for(uint64_t j = 0; j < 3; ++j) {
            result->components.u[i] += metric[i][j] * vector->components.u[j];
        }
    }

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
    _getMetricAt(manifold, manifold->metricField, coords, &metric);
    RmcFloat length = 0;
    for(uint64_t i = 0; i < 3; ++i) {
        for(uint64_t j = 0; j < 3; ++j) {
            length += metric[i][j] * vector->components.u[i] * vector->components.u[j];
        }
    }
    *result = sqrt(length);
    return RMC_SUCCESS;
}