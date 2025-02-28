#include <cmath>
#include <stdlib.h>
#include <raylib.h>
#include <unistd.h>
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <mpi.h>

#include "common.h"
#include "Particle.h"
#include "equations.h"
#include "io.h"

void
compute_acceleration_for_bucket(Particle_Bucket* bucket, Particle* part_i, double smooth)
{
    assert(part_i);
    assert(!part_i->depend);

    if (bucket) {
        Particle* part_j = *bucket;
        while (part_j) {
            if (part_j != part_i)
                calculate_acceleration(part_i, part_j, smooth);

            assert(part_j->next != part_j);
            part_j = part_j->next;
        }
    }
}

void
one_step(Particles* particles, Particle* part_arr, Simulator_Params* params, Extent partition_extent)
{
    for (int x = 0; x < params->n_cells_x; x++) {
        for (int y = 0; y < params->n_cells_y; y++) {
            Particle_Bucket* this_bucket = get_bucket(params, particles, x, y);
            assert(this_bucket);

            Particle* part_i;

            if ((part_i = *this_bucket)) {
                while (part_i) {
                    part_i->acceleration.x = 0; part_i->acceleration.y = 0;

                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            if (!part_i->depend)
                                compute_acceleration_for_bucket(get_bucket(params, particles, x+i-1, y+i-1), part_i, params->smoothing);
                    part_i = part_i->next;
                }
            }
        }
    }

    for (int i = 0; i < params->n_particles; i++) {
        if (!part_arr[i].invalid) {
            part_arr[i].velocity.x += params->time_step * part_arr[i].acceleration.x;
            part_arr[i].velocity.y += params->time_step * part_arr[i].acceleration.y;

            double new_x = part_arr[i].position.x + params->time_step * part_arr[i].velocity.x;
            double new_y = part_arr[i].position.y + params->time_step * part_arr[i].velocity.y;
            if (IN_EXTENT_X(new_x, partition_extent) && IN_EXTENT_Y(new_y, partition_extent))
                move_particle(particles, &part_arr[i], params, partition_extent, new_x, new_y);
            else
                remove_particle(particles, &part_arr[i], params, partition_extent);
        }

    }
}

int
main(int argc, char** argv)
{
    Simulator_Params params = {
        .n_particles = 3,
        .time_step = 0.1,
        .seed = 1000,

        .screen_x = 1920,
        .screen_y = 1080,

        .object_radius = 10,
        .smoothing = 1,

        .n_cells_x = 16,
        .n_cells_y = 9,
        .grid_length = 120,
    };

    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    Extent partition_extent, depend_extent;
    sync_init(&params, &partition_extent, &depend_extent);

    srand(params.seed);
    Particle* particles = new Particle[params.n_particles*4];

    Particles part_hash_map;
    part_hash_map.buckets = new Particle*[params.n_cells_x * params.n_cells_y];

    for (int i = 0; i < params.n_cells_y * params.n_cells_x; i++)
        part_hash_map.buckets[i] = NULL;

    for (int i = 0; i < params.n_particles*4; i++) {
        particles[i].invalid = true;
    }

    for (int i = params.n_particles*rank; i < params.n_particles; i++) {
        particles[i].position.x = static_cast<double>((rand() % static_cast<int>(partition_extent.w)) + partition_extent.x);
        particles[i].position.y = static_cast<double>((rand() % static_cast<int>(partition_extent.h)) + partition_extent.y);
        particles[i].velocity.x = 0;
        particles[i].velocity.y = 0;
        particles[i].prev = NULL;
        particles[i].next = NULL;
        particles[i].invalid = false;
        particles[i].depend = false;

        insert_particle(&part_hash_map, &particles[i], &params, partition_extent);
        particles[i].mass = rand();
    }

    for (int i = 0; i < 10000000; i++) {
        one_step(&part_hash_map, particles, &params, partition_extent);
        sync_all(rank, particles, rank * params.n_particles,
                params.n_particles, &part_hash_map, depend_extent,
                partition_extent, &params);
        printf("%d\n", i);
    }

    MPI_Finalize();

    return 0;
}
