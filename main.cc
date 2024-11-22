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
    // bucket_cache[y][x];
    Particle_Bucket* bucket_cache[3][3];
    for (int x = 0; x < params->n_cells_x; x++) {
        for (int y = 0; y < params->n_cells_y; y++) {
            bucket_cache[0][0] = get_bucket(params, particles, x-1, y-1);
            bucket_cache[0][1] = get_bucket(params, particles, x, y-1);
            bucket_cache[0][2] = get_bucket(params, particles, x+1, y-1);

            bucket_cache[1][0] = get_bucket(params, particles, x-1, y);
            bucket_cache[1][1] = get_bucket(params, particles, x, y);
            bucket_cache[1][2] = get_bucket(params, particles, x+1, y);

            bucket_cache[2][0] = get_bucket(params, particles, x-1, y+1);
            bucket_cache[2][1] = get_bucket(params, particles, x, y+1);
            bucket_cache[2][2] = get_bucket(params, particles, x+1, y+1);

            assert(bucket_cache[1][1]);

            Particle* part_i;

            if ((part_i = *bucket_cache[1][1])) {
                while (part_i) {
                    part_i->acceleration.x = 0; part_i->acceleration.y = 0;

                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            compute_acceleration_for_bucket(bucket_cache[i][j], part_i, params->smoothing);
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

void
print_particle(Particle* particles, int j, int i)
{
    printf("%d: mass = %f, acc = <%f, %f>, vel = <%f, %f>, pos = (%f, %f)\n",
            i, particles[j].mass,
            particles[j].acceleration.x, particles[j].acceleration.y,
            particles[j].velocity.x, particles[j].velocity.y,
            particles[j].position.x, particles[j].position.y);
}

int
main(void)
{
    Simulator_Params params = {
        .n_particles = 3,
        .time_step = 0.1,
        .seed = 1000,

        .screen_x = 1920,
        .screen_y = 1080,

        .object_radius = 20,
        .smoothing = 0.1,

        .render_after_n_frames = 20,

        .n_cells_x = 1,
        .n_cells_y = 1,
        .grid_length = 1000
    };

    srand(params.seed);
    Particle* particles = new Particle[params.n_particles];
    Particles part_hash_map;
    part_hash_map.buckets = new Particle*[params.n_cells_x * params.n_cells_y];
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

    InitWindow(params.screen_x, params.screen_y, "nbody");
    SetTargetFPS(10000);

    auto tm = std::chrono::high_resolution_clock::now();
    auto dur = tm.time_since_epoch();

    while (!WindowShouldClose()) {
//    for (uint64_t i = 0; i < 100000000; i++) {

       BeginDrawing();
       ClearBackground(RAYWHITE);
        for (int j = 0; j < params.n_particles; j++) {
            DrawCircle((int)particles[j].position.x, (int)particles[j].position.y, params.object_radius, DARKBLUE);
        }

        EndDrawing();

        for (int i = 0; i < params.render_after_n_frames; i++)
            one_step(&part_hash_map, particles, &params);
    }

    CloseWindow();
    return 0;
}
