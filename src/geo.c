#include "geo.h"
#include "rmc.h"

RmcError rmcTensorSum(const RmcTensor *a, const RmcTensor *b, RmcTensor *result) {
    if(!a || !b || !result) {
        return RMC_ERROR_GENERIC_NULLPTR;
    }

    if(a->type != b->type) {
        return RMC_ERROR_TENSOR_TYPES_INCOMPATIBLE;
    }

    *result = (RmcTensor){
        .type = a->type,
        .components = {
            .x = a->components.x + b->components.x,
            .y = a->components.y + b->components.y,
            .z = a->components.z + b->components.z
        }
    };

    return RMC_SUCCESS;    
}

RmcError rmcTensorMul(RmcFloat r, const RmcTensor *t, RmcTensor *result) {
    if(!t || !result) {
        return RMC_ERROR_GENERIC_NULLPTR;
    }

    *result = (RmcTensor) {
        .type = t->type,
        .components= {
            .x = r * t->components.x,
            .y = r * t->components.y,
            .z = r * t->components.z
        }
    };

    return RMC_SUCCESS;
}