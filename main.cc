#include <cmath>
#include <stdlib.h>
#include <raylib.h>
#include <unistd.h>
#include <stdio.h>
#include <chrono>
#include <assert.h>
#include <stdint.h>

#include "common.h"
#include "Particle.h"
#include "equations.h"

void
compute_acceleration_for_bucket(Particle_Bucket* bucket, Particle* part_i, double smooth)
{
    assert(part_i);
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
one_step(Particles* particles, Particle* part_arr, Simulator_Params* params)
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
                            compute_acceleration_for_bucket(get_bucket(params, particles, x+i-1, y+i-1), part_i, params->smoothing);
                    part_i = part_i->next;
                }
            }
        }
    }

    for (int i = 0; i < params->n_particles; i++) {
        // if particle is not in a grid ignore it
        int x_index = BUCKET_X(part_arr[i].position.x);
        int y_index = BUCKET_Y(part_arr[i].position.y);

        if (0 <= x_index && x_index < params->n_cells_x && 0 <= y_index && y_index < params->n_cells_y) {
            part_arr[i].velocity.x += params->time_step * part_arr[i].acceleration.x;
            part_arr[i].velocity.y += params->time_step * part_arr[i].acceleration.y;

            double new_x = part_arr[i].position.x + params->time_step * part_arr[i].velocity.x;
            double new_y = part_arr[i].position.y + params->time_step * part_arr[i].velocity.y;
            move_particle(particles, &part_arr[i], params, new_x, new_y);
        }

    }
}

int
main(void)
{
    Simulator_Params params = {
        .n_particles = 2000,
        .time_step = 0.1,
        .seed = 1000,

        .screen_x = 1920,
        .screen_y = 1080,

        .object_radius = 10,
        .smoothing = 1,

        .n_cells_x = 16,
        .n_cells_y = 9,
        .grid_length = 150
    };

    srand(params.seed);
    Particle* particles = new Particle[params.n_particles];

    Particles part_hash_map;
    part_hash_map.buckets = new Particle*[params.n_cells_x * params.n_cells_y];

    InitWindow(params.screen_x, params.screen_y, "nbody");
    SetTargetFPS(10000);

    for (int i = 0; i < params.n_cells_y*params.n_cells_x; i++)
        part_hash_map.buckets[i] = NULL;

    for (int i = 0; i < params.n_particles; i++) {
        particles[i].position.x = static_cast<double>(rand() % (static_cast<int>(params.grid_length) * params.n_cells_x));
        particles[i].position.y = static_cast<double>(rand() % (static_cast<int>(params.grid_length) * params.n_cells_y));
        particles[i].velocity.x = 0;
        particles[i].velocity.y = 0;
        particles[i].prev = NULL;
        particles[i].next = NULL;

        insert_particle(&part_hash_map, &particles[i], &params);
        particles[i].mass = rand();
    }

    for (int i = 0; i < 60000; i++) {
        BeginDrawing();
        ClearBackground(BLACK);

        for (int j = 0; j < params.n_particles; j++)
            DrawCircle((int)particles[j].position.x, (int)particles[j].position.y, params.object_radius, DARKBLUE);

        EndDrawing();

        one_step(&part_hash_map, particles, &params);
    }

    CloseWindow();
    return 0;
}
