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

    RmcTensor a = {
        RMC_TENSOR_TYPE_VECTOR,
        {
            .x = 3,
            .y = 4,
            .z = 0
        }
    };

    RmcFloat result;
    rmcManifoldTensorGetLength(manifold, &a.components, &a, &result);
    printf("Length of vector is %f\n", result);

    rmcManifoldDestroy(manifold);

    return 0;
}