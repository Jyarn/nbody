#pragma once

#include "Particle.h"
#include "common.h"

#define SYNC_TERM (int)(-1)
#define SYNC_TAG 0

// define ordering of doubles when syncing particles between nodes
#define SYNC_POS_X 0
#define SYNC_POS_Y 1
#define SYNC_VEL_X 2
#define SYNC_VEL_Y 3
#define SYNC_ACC_X 4
#define SYNC_ACC_Y 5

void sync_init(int rank, Simulator_Params* params, Extent* depend_extent_ret,
        Extent* recv_extent_ret);

int sync_all(int rank, Particle* part_arr, int particle_partition_start,
        int particle_partiton_len, Particles* hash_map, Extent depend_extent,
        Extent partition_extent, Simulator_Params* params);
