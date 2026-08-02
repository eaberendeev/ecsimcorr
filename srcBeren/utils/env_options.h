#pragma once

#include "util.h"

namespace envOptions {
static inline int validationPeriodicity() {
    const static int res = getenvParsed<int>("VALIDATION_PERIODICITY", 10);
    return res;
}

static inline int timeredSimSteps() {
    const static int res = getenvParsed<int>("TIMERED_SIM_STEPS", 20);
    return res;
}

static inline bool useMixedPrecision() {
    const static bool res = getenvParsed<bool>("USE_MIXED_PRECISION", true);
    return res;
}

}   // namespace envOptions
