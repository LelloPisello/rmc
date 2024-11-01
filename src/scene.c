#include "scene.h"
#include "rmc.h"
#include <stdlib.h>

/*
    SCENE
*/

struct RmcScene_s {
    uint32_t uVolumeCount;
    RmcVolumeInfo *pVolumes;
};

RmcError rmcSceneCreate(const RmcSceneCreateInfo *info, RmcScene *result) {
    if(!info || !result) {
        return RMC_ERROR_GENERIC_NULLPTR;
    }

    if(info->uVolumeCount && !info->pVolumes) {
        return RMC_ERROR_GENERIC_NULLPTR;
    }

#define ALIAS (*result)
    ALIAS = malloc(sizeof(*ALIAS));
    if(!ALIAS) {
        return RMC_ERROR_GENERIC_FAILED_ALLOCATION;
    }
    ALIAS->uVolumeCount = info->uVolumeCount;

    if(info->uVolumeCount) {
        ALIAS->pVolumes = calloc(ALIAS->uVolumeCount, sizeof(RmcVolumeInfo));
        if(!ALIAS->pVolumes) {
            free(ALIAS);
            return RMC_ERROR_GENERIC_FAILED_ALLOCATION;
        }
        for(uint32_t i = 0; i < info->uVolumeCount; ++i) {
            ALIAS->pVolumes[i] = info->pVolumes[i];
        }
    } else {
        ALIAS->pVolumes = NULL;
    }

#undef ALIAS

    return RMC_SUCCESS;
}

void rmcSceneDestroy(RmcScene scene) {
    free(scene->pVolumes);
    free(scene);
}

/*
    RAY
*/