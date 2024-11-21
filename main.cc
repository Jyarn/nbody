#include <cmath>
#include <stdlib.h>
#include <raylib.h>
#include <unistd.h>
#include <stdio.h>
#include <chrono>
#include <string>

#define G 6.67e-11

typedef struct {
    double x, y;
} Vector;

typedef struct {
    Vector velocity, position, acceleration;
    double mass;
} Particle;

typedef struct {
    int n_particles;
    int n_steps;
    double time_step;
    int seed;
    int screen_x, screen_y;
    int obj_x_max, obj_y_max;
    int object_radius;
    double smoothing;
    double clamp_radius;
    int render_after_n_frames;
    double weight_multiplier;
    double discard;
} Simulator_Params;

double
norm(Vector i)
{
    return std::sqrt(i.x*i.x + i.y*i.y);
}

void
one_step(Particle* particles, const Simulator_Params* params)
{
    for (int i = 0; i < params->n_particles; i++) {
        particles[i].acceleration.x = 0;
        particles[i].acceleration.y = 0;
        for (int j = 0; j < params->n_particles; j++) {
            if (i != j) {
                Particle part_i = particles[i];
                Particle part_j = particles[j];

                Vector r;
                r.x = part_j.position.x - part_i.position.x;
                r.y = part_j.position.y - part_i.position.y ;
                if (r.x > params->discard && r.y > params->discard)
                    continue;

                double mass = part_j.mass;
                double r_norm = norm(r) + params->smoothing;
                double r_norm_cubed = r_norm * r_norm * r_norm;

                particles[i].acceleration.x += G * mass * r.x / r_norm_cubed;
                particles[i].acceleration.y += G * mass * r.y / r_norm_cubed;
            }
        }
    }

    for (int i = 0; i < params->n_particles; i++) {
        particles[i].velocity.x += params->time_step * particles[i].acceleration.x;
        particles[i].velocity.y += params->time_step * particles[i].acceleration.y;

        particles[i].position.x += params->time_step * particles[i].velocity.x;
        particles[i].position.y += params->time_step * particles[i].velocity.y;
    }
}

void
print_particle(Particle* particles, int j, int i)
{
    if (norm(particles[j].acceleration) > 0.1)
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
        .obj_x_max = 1920,
        .obj_y_max = 1080,
        .object_radius = 20,
        .smoothing = 10,
        .render_after_n_frames = 20,
        .discard = 1000,
    };

    srand(params.seed);
    Particle* particles = new Particle[params.n_particles];

    for (int i = 0; i < params.n_particles; i++) {
        particles[i].position.x = rand() % params.obj_x_max - (double)params.obj_x_max / 2 + (double)params.screen_x / 2;
        particles[i].position.y = rand() % params.obj_y_max - (double)params.obj_y_max / 2 + (double)params.screen_y / 2;
        particles[i].mass = rand();
    }

    InitWindow(params.screen_x, params.screen_y, "nbody");
    SetTargetFPS(10000);

    auto tm = std::chrono::high_resolution_clock::now();
    auto dur = tm.time_since_epoch();

    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(RAYWHITE);
        for (int j = 0; j < params.n_particles; j++) {
            DrawCircle((int)particles[j].position.x, (int)particles[j].position.y, params.object_radius, DARKBLUE);
        }

        EndDrawing();

        for (int i = 0; i < params.render_after_n_frames; i++)
            one_step(particles, &params);
    }

    CloseWindow();
    return 0;
}
