#include <stdio.h>
#include "manifold.h"
#include "rmc.h"

int main(int argC, char* argV[]) {
    RmcTensor a = {RMC_TENSOR_TYPE_VECTOR, 
    1, 2, 3};
    RmcTensor b = {RMC_TENSOR_TYPE_VECTOR,
    0, 1, 5};

    rmcTensorSum(&a, &b, &b);
    rmcTensorMul(0.5, &b, &b);

    RmcManifold manifold;

    {
        RmcManifoldCreateInfo info = {
            .fBounds = 1.0,
            .uResolution = 32,
            .fpMetricGenerator = (RmcMetricGenerator)1,
        };
        rmcManifoldCreate(&info, &manifold);
    }

    printf("%f %f\n", b.components.x, b.components.y);

    rmcManifoldDestroy(manifold);

    return 0;
}