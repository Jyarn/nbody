#include "common.h"
#include "Particle.h"

#include <cmath>

#define G 6.67e-11

double
norm(Vector i)
{
    return std::sqrt(i.x*i.x + i.y*i.y);
}


void
calculate_acceleration(Particle* part_i, Particle* part_j, double smoothing)
{
    Vector r;
    r.x = part_j->position.x - part_i->position.x;
    r.y = part_j->position.y - part_i->position.y ;

    double mass = part_j->mass;
    double r_norm = norm(r) + smoothing;
    double r_norm_cubed = r_norm * r_norm * r_norm;

    part_i->acceleration.x += G * mass * r.x / r_norm_cubed;
    part_i->acceleration.y += G * mass * r.y / r_norm_cubed;
}
