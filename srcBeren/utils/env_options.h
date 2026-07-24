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

}   // namespace envOptions
