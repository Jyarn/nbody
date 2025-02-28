#pragma once
#include <stdbool.h>

#include "common.h"

// assumes we have a variable of type Simulator_Params* named params in scope
#define BUCKET_X(x) (static_cast<int>(x) / params->grid_length)
#define BUCKET_Y(y) (static_cast<int>(y) / params->grid_length)

typedef struct particle {
    Vector velocity, position, acceleration;
    double mass;
    bool invalid;
    bool depend;
    struct particle* next, *prev;
} Particle;

typedef Particle* Particle_Bucket;

typedef struct {
    Particle_Bucket* buckets;
} Particles;


Particle_Bucket* get_bucket(Simulator_Params* params, Particles* particles,
        int x, int y);

Particle_Bucket* get_particle_bucket(Particles* particles, Particle* part,
        Simulator_Params* params, Extent partition_extent);

int insert_particle(Particles* particles, Particle* part,
        Simulator_Params* params, Extent partition_extent);

void remove_particle(Particles* particles, Particle* part,
        Simulator_Params* params, Extent partition_extent);

void move_particle(Particles* particles, Particle* part,
        Simulator_Params* params, Extent partition_extent,
        double new_x, double new_y);
