#include "Particle.h"
#include "common.h"

#include <assert.h>

#define NULL 0x0

Particle_Bucket*
get_bucket(Simulator_Params* params, Particles* particles, int x, int y)
{
    if ((0 <= x && x < params->n_cells_x)
        && (0<= y && y < params->n_cells_y))
        return &particles->buckets[x + y*params->n_cells_x];
    return NULL;
}

Particle_Bucket*
get_particle_bucket(Particles* particles, Simulator_Params* params, double x, double y)
{
    int x_index = BUCKET_X(x);
    int y_index = BUCKET_Y(y);

    if (0 <= x_index && x_index < params->n_cells_x && 0 <= y_index && y_index < params->n_cells_y)
        return &particles->buckets[x_index + y_index*params->n_cells_x];

    return NULL;
}

int
insert_particle(Particles* particles, Particle* part, Simulator_Params* params)
{
    // check particle is not in another bucket
    assert(particles);
    assert(part && !part->next && !part->prev);

    Particle_Bucket* bucket = get_particle_bucket(particles, params, part->position.x, part->position.y);

    if (!bucket)
        return 1;
    else if (*bucket == NULL) {
        *bucket = part;
    }
    else {
        // check head of bucket is actually the head
        assert(!(*bucket)->prev);
        part->next = *bucket;
        (*bucket)->prev = part;
        *bucket = part;
    }

    return 0;
}

void
remove_particle(Particle* particle)
{
    if (particle->next) {
        // check particle is not head and particle is not the only element in bucket
        assert(particle && particle->next && particle->prev);
        particle->next->prev = particle->prev;
        particle->prev->next = particle->next;

        // check integrity
        assert(particle->next->prev->next == particle->next);
        assert(particle->prev->next->prev == particle->prev);

    } else {
        particle->prev->next = NULL;
    }

   particle->next = NULL;
   particle->prev = NULL;
}

void
move_particle(Particles* particles, Particle* part, Simulator_Params* params, double new_x, double new_y)
{
    int x_index = BUCKET_X(part->position.x);
    int y_index = BUCKET_Y(part->position.y);
    assert(0 <= x_index && y_index < params->n_cells_x);
    assert(0 <= x_index && y_index < params->n_cells_y);

    if ((BUCKET_X(part->position.x) != BUCKET_X(new_x))
        || (BUCKET_Y(part->position.y) != BUCKET_X(new_y)))
    {
        // remove particle from it's bucket
        if (!part->prev) {
            assert((part->next && !part->prev) || !part->next);

            // case: particle is the only one in the bucket or is head
            Particle_Bucket* bucket = get_particle_bucket(particles, params,
                    part->position.x, part->position.y);

            // check particle is head
            assert(bucket && *bucket == part);
            if (part->next)
                part->next->prev = NULL; // <----------- seg fault
            *bucket = part->next;
            part->next = NULL;

        } else {
            // case: particle has other particles in the same bucket
            remove_particle(part);
        }

        // update position
        part->position.x = new_x;
        part->position.y = new_y;

        // insert it into new bucket
        insert_particle(particles, part, params);

    } else {
        // update position
        part->position.x = new_x;
        part->position.y = new_y;
    }
}
