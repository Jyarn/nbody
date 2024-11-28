#pragma once
#include "Particle.h"

void dump_init(int rank);
void dump(Particle* part_arr, int n, Simulator_Params* params);
void dump_clean(void);
