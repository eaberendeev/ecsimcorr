#include "maxwell_operators.h"
#include "timer.h"

void MaxwellOperators::advance_magnetic_field(Field3d& Bn, const Field3d& B, const Field3d& E, const Field3d& En,
                                              double dt) const {
    RECORD_TIMER;
    Bn = B - (0.5 * dt) * (curlE * (E + En));
}
