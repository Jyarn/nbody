#pragma once

#include <assert.h>

typedef struct {
    double x, y;
} Vector;

typedef struct {
    int n_particles;
    double time_step;
    int seed;
    int screen_x, screen_y;
    int object_radius;
    double smoothing;
    int n_cells_x, n_cells_y;
    double grid_length;
    int n_iterations;
    int self_rank;
} Simulator_Params;

typedef struct {
    double w, h, x, y;
} Extent;

#define IN_EXTENT_X(xx, extent) (extent.x <= xx && xx < extent.x + extent.w)
#define IN_EXTENT_Y(yy, extent) (extent.y <= yy && yy < extent.y + extent.h)
