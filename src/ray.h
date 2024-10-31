#ifndef RMC_RAY_h_
#define RMC_RAY_h_

#include "geo.h"
#include "manifold.h"
#include "rmc.h"

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
} RmcRayMarchParams;

typedef struct {
    RmcRayHitType type;
    RmcFloat fDistance;
    RmcColor color;
    uint64_t uStepCount;
} RmcRayHitInfo;

RmcError rmcRayMarch(RmcManifold manifold, const RmcRayMarchParams* info, RmcRay* ray, RmcRayHitInfo* result);

#endif