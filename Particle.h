#pragma once
#include "common.h"

// assumes we have a variable of type Simulator_Params* named params in scope
#define BUCKET_X(x) (static_cast<int>(x) / params->grid_length)
#define BUCKET_Y(y) (static_cast<int>(y) / params->grid_length)

typedef struct particle {
    Vector velocity, position, acceleration;
    double mass;
    struct particle* next, *prev;
} Particle;

typedef Particle* Particle_Bucket;

typedef struct {
    Particle_Bucket* buckets;
} Particles;


Particle_Bucket* get_bucket(Simulator_Params* params, Particles* particles,
        int x, int y);

Particle_Bucket* get_particle_bucket(Particles* particles,
        Simulator_Params* params, double x, double y);

int insert_particle(Particles* particles, Particle* part,
        Simulator_Params* params);

void remove_particle(Particle* particle);

void move_particle(Particles* particles, Particle* part,
        Simulator_Params* params, double new_x, double new_y);
