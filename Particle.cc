#include "Particle.h"
#include "common.h"

#include <stdlib.h>
#include <string.h>

#define ARR_THRESHOLD 2

Particle_Bucket*
get_bucket(Simulator_Params* params, Particles* particles, int x, int y)
{
    if ((0 <= x && x < params->n_cells_x)
        && (0 <= y && y < params->n_cells_y))
        return &particles->buckets[x + y*params->n_cells_x];
    return NULL;
}

Particle_Bucket*
get_particle_bucket(Particles* particles, Particle* part,
        Simulator_Params* params, Extent partition_extent)
{
    int x_index = BUCKET_X(part->position.x);
    int y_index = BUCKET_Y(part->position.y);
    assert(0 <= x_index && x_index < params->n_cells_x);
    assert(0 <= y_index && y_index < params->n_cells_y);

    return &particles->buckets[x_index + y_index*params->n_cells_x];
}

int
insert_particle(Particles* particles, Particle* part, Simulator_Params* params,
        Extent partition_extent)
{
    // check particle is not in another bucket
    assert(particles);
    assert(part->invalid);
    assert(part->part_idx == -1);

    // bypass assert
    Particle_Bucket* bucket = get_particle_bucket(particles, part, params,
            partition_extent);

    assert(bucket);
    assert(bucket->size <= bucket->length);

    // resize bucket
    if (bucket->size == bucket->length) {
        assert(bucket->length != 0);
        Particle** new_bucket = new Particle*[bucket->length *= ARR_THRESHOLD];
        memcpy(new_bucket, bucket->arr, sizeof(Particle*) * (bucket->size));
        bucket->arr = new_bucket;
    }

    assert(bucket->size < bucket->length);
    bucket->arr[bucket->size] = part;
    part->part_idx = bucket->size++;
    part->invalid = false;
    return 0;
}

void
remove_particle(Particles* particles, Particle* part, Simulator_Params* params,
        Extent partition_extent)
{
    Particle_Bucket* bucket = get_particle_bucket(particles, part, params,
            partition_extent);
    int part_idx = part->part_idx;

    assert(!part->invalid);
    assert(part == bucket->arr[part_idx]);
    assert(bucket->size <= bucket->length);

    // shift left
    for (int i = part_idx; i < bucket->size - 1; i++) {
        assert(i > bucket->length || bucket->arr[i+1]->part_idx == i+1);
        bucket->arr[i] = bucket->arr[i+1];
        bucket->arr[i]->part_idx = i;
    }

    assert(bucket->size <= 1 || part_idx == bucket->size -1 || bucket->arr[bucket->size - 2] == bucket->arr[bucket->size-1]);

    // trim arr
    if (bucket->size < bucket->length / ARR_THRESHOLD) {
        assert(bucket->length / ARR_THRESHOLD != 0);
        Particle** new_bucket = new Particle*[bucket->length /= ARR_THRESHOLD];
        memcpy(new_bucket, bucket->arr, sizeof(Particle*) * (--bucket->size));

    } else {
        bucket->size--;
    }

    part->part_idx = -1;
    part->invalid = true;
}

void
move_particle(Particles* particles, Particle* part, Simulator_Params* params,
        Extent partition_extent, double new_x, double new_y)
{
    int x_index = BUCKET_X(part->position.x);
    int y_index = BUCKET_Y(part->position.y);

    assert(0 <= x_index && x_index < params->n_cells_x);
    assert(0 <= y_index && y_index < params->n_cells_y);
    assert(IN_EXTENT_X(part->position.x, partition_extent));
    assert(IN_EXTENT_Y(part->position.y, partition_extent));
    assert(IN_EXTENT_X(new_x, partition_extent));
    assert(IN_EXTENT_Y(new_y, partition_extent));

    assert(!part->invalid);
    assert(part->part_idx != -1);

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

void
init_particles(Particle* part_arr, Particles* part_hash_map, Extent partition_extent, Simulator_Params* params)
{
    int part_start = params->self_rank * params->n_particles;
    int part_end = part_start + params->n_particles;

    srand(params->seed * (params->self_rank + 1));
    part_hash_map->buckets = new Particle_Bucket[params->n_cells_x * params->n_cells_y];

    for (int i = 0; i < params->n_cells_y * params->n_cells_x; i++) {
        part_hash_map->buckets[i].size = 0;
        part_hash_map->buckets[i].length = 1;
        part_hash_map->buckets[i].arr = new Particle*[1];
    }

    for (int i = 0; i < params->n_particles * 4; i++) {
        if (part_start <= i && i < part_end) {
            part_arr[i].mass = rand();
            part_arr[i].position.x = static_cast<double>((rand() % static_cast<int>(partition_extent.w)) + partition_extent.x);
            part_arr[i].position.y = static_cast<double>((rand() % static_cast<int>(partition_extent.h)) + partition_extent.y);
            part_arr[i].velocity.x = 0;
            part_arr[i].velocity.y = 0;
            part_arr[i].invalid = true;
            part_arr[i].depend = false;
            part_arr[i].part_idx = -1;

            insert_particle(part_hash_map, &part_arr[i], params,
                    partition_extent);
            part_arr[i].mass = rand();

        } else {
            part_arr[i].position.x = -1;
            part_arr[i].position.y = -1;
            part_arr[i].velocity.x = 0;
            part_arr[i].velocity.y = 0;
            part_arr[i].invalid = true;
            part_arr[i].depend = false;
            part_arr[i].part_idx = -1;
            part_arr[i].mass = 0;
        }
    }
}
