#ifndef RMC_RAY_h_
#define RMC_RAY_h_

#include "geo.h"
#include "manifold.h"
#include "rmc.h"

//VOLUMES

typedef struct RmcVolumeInfo_s RmcVolumeInfo;

//hit is a boolean value, extra data can be freely interpreted
typedef void (*RmcVolumeFunction)(void* pData, const RmcCoordinates* coords, RmcBool *hit);

//i am aware this sucks thanks
struct RmcVolumeInfo_s {    
    RmcVolumeFunction pfVolume;
    void* pData;
};

//SCENE

typedef struct {
    uint32_t uVolumeCount;
    RmcVolumeInfo *pVolumes;
} RmcSceneCreateInfo;

//i chose to make the scene an opaque handle instead of a simple list
//to make it possible to restructure easily later 
RmcError rmcSceneCreate(const RmcSceneCreateInfo* info, RmcScene *result);
void rmcSceneDestroy(RmcScene);

//RAYS

typedef enum {
    RMC_RAY_HIT_OBJECT,
    RMC_RAY_HIT_BOUNDS,
    //when exceeding limits
    RMC_RAY_HIT_MISSED,
} RmcRayHitType;

typedef union {
    RmcFloat u[3];
    struct {
        RmcFloat r, g, b;
    };
} RmcColor;

typedef struct {
    RmcTensor direction;
    RmcCoordinates position;
} RmcRay;

typedef struct {
    uint64_t uMaxSteps;
    //ray.dir * stepSize
    RmcFloat fStepSize;
    //set to < 0 for infinite
    RmcFloat fMaxDistance;
} RmcRayMarchInfo;

typedef struct {
    RmcRayHitType type;
    RmcFloat fDistance;
    RmcColor color;
    uint64_t uStepCount;
} RmcRayHitInfo;

RmcError rmcRayMarch(
    RmcManifold manifold, 
    RmcScene scene, 
    const RmcRayMarchInfo* info, 
    RmcRay* ray, 
    RmcRayHitInfo* result
);

#endif