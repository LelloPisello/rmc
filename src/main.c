#include <stdio.h>
#include "geo.h"
#include "manifold.h"
#include "rmc.h"
#include "scene.h"

//kroenecker delta
void metrFill(const RmcCoordinates* coords, RmcMetricOutput* output) {
    for(uint64_t i = 0; i < 3; ++i) {
        for(uint64_t j = 0; j < 3; ++j) {
            (*output)[i][j] = (RmcFloat)(i == j);
        }
    }
}

void sphere(void *pData, const RmcCoordinates* coords, RmcBool* hit) {
    *hit = 0;
    if(coords->x * coords->x +
    coords->y * coords->y +
    coords->z * coords->z < 1.0) {
        *hit = 1;
    }
}

int main(int argC, char* argV[]) {

    RmcManifold manifold;
    RmcScene scene;

    {
        RmcManifoldCreateInfo info = {
            .fBounds = 1.0,
            .uMetricResolution = 32,
            .fpMetricGenerator = metrFill,
        };
        rmcManifoldCreate(&info, &manifold);
    }

    {
        RmcVolumeInfo volume = {
            .pData = NULL,
            .pfVolume = sphere,
        };
        RmcSceneCreateInfo info = {
            .uVolumeCount = 1,
            .pVolumes = &volume
        };
        rmcSceneCreate(&info, &scene);
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

    printf("Raising index of tensor\n");
    rmcManifoldTensorRaiseIndex(manifold, &coords, &a, &b);

    printf("Getting tensor length\n");
    RmcFloat result;
    rmcManifoldTensorGetLength(manifold, &coords, &b, &result);
    printf("Length of vector is %f\n", result);

    rmcManifoldDestroy(manifold);
    rmcSceneDestroy(scene);

    return 0;
}