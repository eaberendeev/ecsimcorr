#ifndef PARTICLES_DISTRIBUTION_COLLECTION_H_
#define PARTICLES_DISTRIBUTION_COLLECTION_H_
#include <vector>

#include "Particle.h"
#include "Vec.h"
#include "random_generator.h"
#include "service.h"
#include "util.h"

void distribute_uniform_cilinderZ(std::vector<Particle>& particles,
                                  const int count, const double3& c,
                                  const double r0, const double z0,
                                  ThreadRandomGenerator& randGenSpace);
void distribute_uniform_rectangle(std::vector<Particle>& particles,
                                  const int count, const double3& c,
                                  const double3& length,
                                  ThreadRandomGenerator& randGenSpace);
void distribute_pulse_gauss(std::vector<Particle>& particles,
                            const double3& sigma,
                            ThreadRandomGenerator& randGenPulse);

#endif
