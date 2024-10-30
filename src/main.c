#include <stdio.h>
#include "geo.h"
#include "manifold.h"
#include "rmc.h"

void metrFill(const RmcCoordinates* coords, RmcMetricOutput* output) {
    for(uint32_t i = 0; i < 3; ++i) {
        for(uint32_t j = 0; j < 3; ++j) {
            (*output)[i][j] = (RmcFloat)(i == j);
        }
    }
}

int main(int argC, char* argV[]) {

    RmcManifold manifold;

    {
        RmcManifoldCreateInfo info = {
            .fBounds = 1.0,
            .uResolution = 32,
            .fpMetricGenerator = metrFill,
        };
        rmcManifoldCreate(&info, &manifold);
    }
    RmcTensor b;
    RmcTensor a = {
        RMC_TENSOR_TYPE_COVECTOR,
        {
            .x = 3,
            .y = 4,
            .z = 0
        }
    };

    RmcCoordinates coords = {
        .x = 0,
        .y = 0,
        .z = 0
    };

    rmcManifoldTensorRaiseIndex(manifold, &coords, &a, &b);

    RmcFloat result;
    rmcManifoldTensorGetLength(manifold, &coords, &b, &result);
    printf("Length of vector is %f\n", result);

    rmcManifoldDestroy(manifold);

    return 0;
}