#pragma once

#include "Particle.h"
#include "common.h"

#define SYNC_TERM (int)(-1)
#define SYNC_TAG 0

#define SYNC_BUFFER_SIZE 7

// define ordering of doubles when syncing particles between nodes
#define SYNC_POS_X 0
#define SYNC_POS_Y 1
#define SYNC_VEL_X 2
#define SYNC_VEL_Y 3
#define SYNC_ACC_X 4
#define SYNC_ACC_Y 5
#define SYNC_MASS 6


int sync_all(Particle* part_arr, Particles* hash_map,
        Simulator_Params* params);

void sync_init(Simulator_Params* params, Extent* depend_extent_ret,
        Extent* partition_extent_ret);
