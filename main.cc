#include <cmath>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <mpi.h>
#include <stdio.h>
#include <omp.h>

#include "common.h"
#include "Particle.h"
#include "io.h"

#define G 6.67e-11

void
compute_acceleration_for_bucket(Particle_Bucket* bucket, Particle* part_i, double smooth)
{
    assert(part_i);
    assert(!part_i->depend);

    if (bucket) {
        double acc_x = 0; double acc_y = 0;
        double part_i_x = part_i->position.x;
        double part_i_y = part_i->position.y;

        #pragma omp parallel for default(none) shared(bucket) firstprivate(part_i_x, part_i_y, smooth) reduction(+:acc_x, acc_y)
        for (int j = 0; j < bucket->size; j++) {
            Particle* part_j = bucket->arr[j];
            Vector r;
            r.x = part_j->position.x - part_i_x;
            r.y = part_j->position.y - part_i_y;

            assert(r.x <= 240);
            assert(r.y <= 240);

            double mass = part_j->mass;
            double r_norm = std::sqrt(r.x*r.x + r.y*r.y) + smooth;
            double r_norm_cubed = r_norm * r_norm * r_norm;

            acc_x += G * mass * r.x / r_norm_cubed;
            acc_y += G * mass * r.y / r_norm_cubed;
        }

        part_i->acceleration.x += acc_x;
        part_i->acceleration.y += acc_y;
    }
}

void
one_step(Particles* particles, Particle* part_arr, Simulator_Params* params, Extent partition_extent)
{
    #pragma omp parallel for default(none) shared(particles) firstprivate(partition_extent, params) collapse(2)
    for (int x = 0; x < params->n_cells_x; x++) {
        for (int y = 0; y < params->n_cells_y; y++) {
            Particle_Bucket* this_bucket = get_bucket(params, particles,
                    x, y);

            assert(this_bucket);

            for (int i = 0; i < this_bucket->size; i++) {
                Particle* part_i = this_bucket->arr[i];
                assert(!part_i->invalid);

                part_i->acceleration.x = 0; part_i->acceleration.y = 0;

                for (int j = 0; j < 3; j++)
                    for (int k = 0; k < 3; k++)
                        if (!part_i->depend)
                            compute_acceleration_for_bucket(
                                    get_bucket(
                                        params, particles, x+k-1, y+j-1
                                        ),
                                    part_i, params->smoothing);
            }
        }
    }

    for (int i = 0; i < params->n_particles; i++) {
        if (!part_arr[i].invalid || part_arr[i].depend) {
            part_arr[i].velocity.x += params->time_step * part_arr[i].acceleration.x;
            part_arr[i].velocity.y += params->time_step * part_arr[i].acceleration.y;

            double new_x = part_arr[i].position.x + params->time_step * part_arr[i].velocity.x;
            double new_y = part_arr[i].position.y + params->time_step * part_arr[i].velocity.y;

            if (IN_EXTENT_X(new_x, partition_extent)
                    && IN_EXTENT_Y(new_y, partition_extent))
            {
                move_particle(particles, &part_arr[i], params,
                        partition_extent, new_x, new_y);

            } else {
                remove_particle(particles, &part_arr[i], params,
                        partition_extent);
            }
        }
    }
}

int
main(int argc, char** argv)
{
    Simulator_Params params;

    params.n_particles = 2000;
    params.time_step = 0.1;
    params.seed = 1000;
    params.screen_x = 1920;
    params.screen_y = 1080;
    params.object_radius = 10;
    params.smoothing = 1;
    params.n_cells_x = 64;
    params.n_cells_y = 36;
    params.grid_length = 120;
    params.n_iterations = 100'000;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &params.self_rank);

    Extent partition_extent, depend_extent;
    sync_init(&params, &depend_extent, &partition_extent);

    Particle* particles = new Particle[params.n_particles*4];

    Particles part_hash_map;

    init_particles(particles, &part_hash_map, partition_extent, &params);


    if (params.self_rank == 0)
        printf("running with %d thread(s)\n", omp_get_max_threads());

    /*
     * I discussed in tutorial with John creating a function that runs this loop
     * to reduce the amount of function calls. I decided not to do this because
     * it ruins readability
     */
    for (int i = 0; i < params.n_iterations; i++) {
        one_step(&part_hash_map, particles, &params, partition_extent);
        sync_all(particles, &part_hash_map, &params);
    }

    MPI_Finalize();
    return 0;
}
