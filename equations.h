#pragma once
#include "common.h"
#include "Particle.h"

double norm(Vector i);
void calculate_acceleration(Particle* part_i, Particle* part_j,
        double smoothing);
