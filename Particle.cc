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
get_particle_bucket(Particles* particles, Particle* part, Simulator_Params* params, Extent partition_extent)
{
    assert(!part->invalid);

    if ((partition_extent.x <= part->position.x && part->position.x < partition_extent.x + partition_extent.w)
        && (partition_extent.y <= part->position.y && part->position.y < partition_extent.y + partition_extent.h))
    {
        int x_index = BUCKET_X(part->position.x);
        int y_index = BUCKET_Y(part->position.y);
        assert(0 <= x_index && x_index < params->n_cells_x && 0 <= y_index && y_index < params->n_cells_y);

        return &particles->buckets[x_index + y_index*params->n_cells_x];
    }

    return NULL;
}

int
insert_particle(Particles* particles, Particle* part, Simulator_Params* params, Extent partition_extent)
{
    // check particle is not in another bucket
    assert(particles);
    assert(part->invalid);
    assert(part && !part->next && !part->prev);

    // bypass assert
    part->invalid = false;
    Particle_Bucket* bucket = get_particle_bucket(particles, part, params, partition_extent);

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

    part->invalid = false;
    return 0;
}

/*
void
remove_particle(Particle* particle)
{
    assert(!particle->invalid);

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
   particle->invalid = true;
}
*/

void
remove_particle(Particles* particles, Particle* part, Simulator_Params* params,
        Extent partition_extent)
{
    // remove particle from it's bucket
    if (!part->prev) {
        assert((part->next && !part->prev) || !part->next);

        // case: particle is the only one in the bucket or is head
        Particle_Bucket* bucket = get_particle_bucket(particles, part, params,
                partition_extent);

        // check particle is head
        assert(bucket && *bucket == part);
        if (part->next)
            part->next->prev = NULL; // <----------- seg fault
        *bucket = part->next;
        part->next = NULL;
        part->invalid = true;

    } else {
        // case: particle has other particles in the same bucket
        assert(!part->invalid);

        if (part->next) {
            // check particle is not head and particle is not the only element in bucket
            assert(part && part->next && part->prev);
            part->next->prev = part->prev;
            part->prev->next = part->next;

            // check integrity
            assert(part->next->prev->next == part->next);
            assert(part->prev->next->prev == part->prev);

        } else {
            part->prev->next = NULL;
        }

        part->next = NULL;
        part->prev = NULL;
        part->invalid = true;
    }

}

void
move_particle(Particles* particles, Particle* part, Simulator_Params* params,
        Extent partition_extent, double new_x, double new_y)
{
    int x_index = BUCKET_X(part->position.x);
    int y_index = BUCKET_Y(part->position.y);

    assert(IN_EXTENT_X(part->position.x, partition_extent));
    assert(IN_EXTENT_Y(part->position.y, partition_extent));
    assert(IN_EXTENT_X(new_x, partition_extent));
    assert(IN_EXTENT_Y(new_y, partition_extent));

    assert(!part->invalid);

    if ((BUCKET_X(part->position.x) != BUCKET_X(new_x))
        || (BUCKET_Y(part->position.y) != BUCKET_X(new_y)))
    {
        remove_particle(particles, part, params, partition_extent);

        // update position
        part->position.x = new_x;
        part->position.y = new_y;

        // insert it into new bucket
        insert_particle(particles, part, params, partition_extent);

    } else {
        // update position
        part->position.x = new_x;
        part->position.y = new_y;
    }
}
