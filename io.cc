#include "io.h"

#include <assert.h>
#include <mpi.h>
#include <stdio.h>

void
send_particle(int index, Particle part, int recv_rank)
{
    MPI_Ssend(&index, 1, MPI_INT, recv_rank, SYNC_TAG, MPI_COMM_WORLD);

    double dbl_buf[6];
    dbl_buf[SYNC_POS_X] = part.position.x;
    dbl_buf[SYNC_POS_Y] = part.position.y;
    dbl_buf[SYNC_VEL_X] = part.velocity.x;
    dbl_buf[SYNC_VEL_Y] = part.velocity.y;
    dbl_buf[SYNC_ACC_X] = part.acceleration.x;
    dbl_buf[SYNC_ACC_Y] = part.acceleration.y;

    MPI_Ssend(dbl_buf, 6, MPI_DOUBLE, recv_rank, SYNC_TAG, MPI_COMM_WORLD);
}

void
sync_terminate_connection(int recv_rank)
{
    int term = SYNC_TERM;
    MPI_Ssend(&term, 1, MPI_INT, recv_rank, SYNC_TAG, MPI_COMM_WORLD);
}

int
sync_send(Particle* part_arr, int n, int recv_rank,
        Extent depend_extent, Extent recv_extent)
{
    Particle part_i;
    Vector part_pos;
    int n_particles_lost = 0;

    for (int i = 0; i < n; i++) {
        part_i = part_arr[i];
        part_pos = part_i.position;

        // skip all valid and dependency particles
        if (!part_i.invalid || part_i.depend) {
            // assert to make sure we aren't skipping particles just because they have certain flags enabled
            assert(!(IN_EXTENT_X(part_pos.x, recv_extent) && IN_EXTENT_Y(part_pos.y, recv_extent)));
            assert((part_i.depend && (IN_EXTENT_X(part_pos.x, depend_extent) && IN_EXTENT_Y(part_pos.y, depend_extent))) || !part_i.depend);
            continue;
        }

        if (IN_EXTENT_X(part_pos.x, recv_extent) && IN_EXTENT_Y(part_pos.y, recv_extent))
        {
            assert(part_i.invalid);
            assert(!part_i.depend);
            assert(!(IN_EXTENT_X(part_pos.x, depend_extent) && IN_EXTENT_Y(part_pos.y, depend_extent)));

            send_particle(i, part_i, recv_rank);
            n_particles_lost++;
        }
    }

    sync_terminate_connection(recv_rank);
    return n_particles_lost;
}

bool
sync_receive_particle(Particle* part_arr, int sender_rank, int* index_ret, double* dbl_ret)
{
    MPI_Status status;
    MPI_Recv(index_ret, 1, MPI_INT, sender_rank, SYNC_TAG, MPI_COMM_WORLD, &status);

    if (*index_ret == SYNC_TERM)
        return false; // return false to terminate loop

    assert(status.MPI_ERROR == MPI_SUCCESS);
    assert(status.MPI_TAG == SYNC_TAG);

    MPI_Recv(dbl_ret, 6, MPI_DOUBLE, sender_rank, SYNC_TAG, MPI_COMM_WORLD, &status);
    assert(status.MPI_ERROR == MPI_SUCCESS);
    assert(status.MPI_TAG == SYNC_TAG);

    return true;
}

int
sync_receive(Particle* part_arr, int particle_partition_start,
        int particle_partition_len, Particles* particle_hash_map,
        int sender_rank, Extent depend_extent, Extent recv_extent,
        Simulator_Params* params)
{
    int part_index = 0;
    double dbl_buf[6] = { 0, 0, 0, 0, 0, 0 };
    int n_particles_gained = 0;

    while (true) {
        if (!sync_receive_particle(part_arr, sender_rank, &part_index, dbl_buf))
            break;

        if (IN_EXTENT_X(dbl_buf[SYNC_POS_X], recv_extent) && IN_EXTENT_Y(dbl_buf[SYNC_POS_Y], recv_extent))
        {
            // assert particle is a valid index in part_arr and part_arr[part_index]
            // is invalid
            assert(particle_partition_start <= part_index);
            assert(part_index < particle_partition_start + particle_partition_len);
            assert(part_arr[part_index].invalid);

            // assert particle is not a dependency of the receiver
            assert(!(IN_EXTENT_X(dbl_buf[SYNC_POS_X], depend_extent) && IN_EXTENT_Y(dbl_buf[SYNC_POS_Y], depend_extent)));

            n_particles_gained++;

            part_arr[part_index].position.x = dbl_buf[SYNC_POS_X];
            part_arr[part_index].position.y = dbl_buf[SYNC_POS_Y];
            part_arr[part_index].velocity.x = dbl_buf[SYNC_VEL_X];
            part_arr[part_index].velocity.y = dbl_buf[SYNC_VEL_Y];
            part_arr[part_index].acceleration.x = dbl_buf[SYNC_ACC_X];
            part_arr[part_index].acceleration.y = dbl_buf[SYNC_ACC_Y];

            insert_particle(particle_hash_map, &part_arr[part_index], params, recv_extent);

        // if particle is in the dependency extent, sync without deleting the particle
        } else if ((depend_extent.x <= dbl_buf[0] && dbl_buf[0] < depend_extent.w + depend_extent.x)
            && (depend_extent.y <= dbl_buf[1] && dbl_buf[1] < depend_extent.h + depend_extent.y))
        {
            // assert that particle is not in recv_extent
            assert(recv_extent.x <= dbl_buf[0] && dbl_buf[0] < recv_extent.w + depend_extent.x);
            assert(recv_extent.y <= dbl_buf[1] && dbl_buf[1] < recv_extent.h + depend_extent.y);

            // assert particle is a valid index in part_arr and part_arr[part_index]
            assert(particle_partition_start <= part_index);
            assert(part_index < particle_partition_start + particle_partition_len);

            // assert(part.invalid => part.depend);
            assert((!part_arr[part_index].invalid && part_arr[part_index].depend) || part_arr[part_index].invalid);

            // mark particle as a dependency and not to calculate acceleration for it during one_step
            part_arr[part_index].depend = true;
            part_arr[part_index].position.x = dbl_buf[SYNC_POS_X];
            part_arr[part_index].position.y = dbl_buf[SYNC_POS_Y];
            part_arr[part_index].velocity.x = dbl_buf[SYNC_VEL_X];
            part_arr[part_index].velocity.y = dbl_buf[SYNC_VEL_Y];
            part_arr[part_index].acceleration.x = dbl_buf[SYNC_ACC_X];
            part_arr[part_index].acceleration.y = dbl_buf[SYNC_ACC_Y];

            if (part_arr[part_index].invalid)
                insert_particle(particle_hash_map, &part_arr[part_index], params, depend_extent);
            else
                move_particle(particle_hash_map, &part_arr[part_index], params, depend_extent, dbl_buf[SYNC_POS_X], dbl_buf[SYNC_POS_Y]);

        } else {
            // equivalent to assert(false)
            assert("Error: node received particle that doesn't belong to any partition" == NULL);
            continue;
        }

    }

    return n_particles_gained;
}

int
sync_all(int rank, Particle* part_arr, int particle_partition_start,
        int particle_partiton_len, Particles* hash_map, Extent depend_extent,
        Extent partition_extent, Simulator_Params* params)
{
    int delta_particles = 0;

    switch (rank) {
        case 0:
            // sync(0, 1)
            delta_particles -= sync_send(part_arr, params->n_particles, 1,
                    depend_extent, partition_extent);
            delta_particles += sync_receive(part_arr, particle_partition_start,
                    particle_partiton_len, hash_map, 1, depend_extent,
                    partition_extent, params);

            // sync(0, 2)
            delta_particles -= sync_send(part_arr, params->n_particles, 2,
                    depend_extent, partition_extent);
            delta_particles += sync_receive(part_arr, particle_partition_start,
                    particle_partiton_len, hash_map, 2, depend_extent,
                    partition_extent, params);

            // sync(0, 3)
            delta_particles -= sync_send(part_arr, params->n_particles, 3,
                    depend_extent, partition_extent);
            delta_particles += sync_receive(part_arr, particle_partition_start,
                    particle_partiton_len, hash_map, 3, depend_extent,
                    partition_extent, params);
           break;
        case 1:
           // sync(1, 0)
           delta_particles += sync_receive(part_arr, particle_partition_start,
                   particle_partiton_len, hash_map, 0, depend_extent,
                   partition_extent, params);
            delta_particles -= sync_send(part_arr, params->n_particles, 0,
                    depend_extent, partition_extent);

            // sync(1, 3)
           delta_particles += sync_receive(part_arr, particle_partition_start,
                   particle_partiton_len, hash_map, 0, depend_extent,
                   partition_extent, params);
            delta_particles -= sync_send(part_arr, params->n_particles, 0,
                    depend_extent, partition_extent);

            // sync(1, 2)
           delta_particles += sync_receive(part_arr, particle_partition_start,
                   particle_partiton_len, hash_map, 2, depend_extent,
                   partition_extent, params);
            delta_particles -= sync_send(part_arr, params->n_particles, 2,
                    depend_extent, partition_extent);
            break;

        case 2:
            // sync(2, 3)
            delta_particles -= sync_send(part_arr, params->n_particles, 3,
                    depend_extent, partition_extent);
            delta_particles += sync_receive(part_arr, particle_partition_start,
                    particle_partiton_len, hash_map, 3, depend_extent,
                    partition_extent, params);

            // sync(2, 0)
           delta_particles += sync_receive(part_arr, particle_partition_start,
                   particle_partiton_len, hash_map, 0, depend_extent,
                   partition_extent, params);
            delta_particles -= sync_send(part_arr, params->n_particles, 0,
                    depend_extent, partition_extent);

            // sync(2, 1)
            delta_particles -= sync_send(part_arr, params->n_particles, 1,
                    depend_extent, partition_extent);
            delta_particles += sync_receive(part_arr, particle_partition_start,
                    particle_partiton_len, hash_map, 1, depend_extent,
                    partition_extent, params);
           break;

        case 3:
            // sync(2, 3)
            delta_particles += sync_receive(part_arr, particle_partition_start,
                    particle_partiton_len, hash_map, 2, depend_extent,
                    partition_extent, params);
            delta_particles -= sync_send(part_arr, params->n_particles, 2,
                    depend_extent, partition_extent);

            // sync(3, 1)
            delta_particles -= sync_send(part_arr, params->n_particles, 1,
                    depend_extent, partition_extent);
           delta_particles += sync_receive(part_arr, particle_partition_start,
                   particle_partiton_len, hash_map, 1, depend_extent,
                   partition_extent, params);

            // sync(3, 0)
            delta_particles += sync_receive(part_arr, particle_partition_start,
                    particle_partiton_len, hash_map, 0, depend_extent,
                    partition_extent, params);
            delta_particles -= sync_send(part_arr, params->n_particles, 0,
                    depend_extent, partition_extent);
           break;
    }

    return delta_particles;
}

void
sync_init(Simulator_Params* params, Extent* depend_extent_ret, Extent* recv_extent_ret)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    switch (rank) {
        case 0:
            recv_extent_ret->x = 0.0;
            recv_extent_ret->y = 0.0;
            recv_extent_ret->w = params->n_cells_x * params->grid_length / 2;
            recv_extent_ret->h = params->n_cells_y * params->grid_length / 2;
            *depend_extent_ret = *recv_extent_ret;
            depend_extent_ret->w += params->grid_length;
            depend_extent_ret->h += params->grid_length;
            break;

        case 1:
            recv_extent_ret->x = 0.0;
            recv_extent_ret->y = params->n_cells_y * params->grid_length / 2;
            recv_extent_ret->w = params->n_cells_x * params->grid_length / 2;
            recv_extent_ret->h = params->n_cells_y * params->grid_length / 2;
            *depend_extent_ret = *recv_extent_ret;
            depend_extent_ret->y -= params->grid_length;
            depend_extent_ret->h += params->grid_length;
            depend_extent_ret->w += params->grid_length;
            break;

        case 2:
            recv_extent_ret->x = params->n_cells_x * params->grid_length / 2;
            recv_extent_ret->y = 0.0;
            recv_extent_ret->w = params->n_cells_x * params->grid_length / 2;
            recv_extent_ret->h = params->n_cells_y * params->grid_length / 2;
            *depend_extent_ret = *recv_extent_ret;
            depend_extent_ret->x -= params->grid_length;
            depend_extent_ret->w += params->grid_length;
            depend_extent_ret->h += params->grid_length;
            break;

        case 3:
            recv_extent_ret->x = params->n_cells_x * params->grid_length / 2;
            recv_extent_ret->y = params->n_cells_y * params->grid_length / 2;
            recv_extent_ret->w = params->n_cells_x * params->grid_length / 2;
            recv_extent_ret->h = params->n_cells_y * params->grid_length / 2;
            *depend_extent_ret = *recv_extent_ret;
            depend_extent_ret->x -= params->grid_length;
            depend_extent_ret->w += params->grid_length;
            depend_extent_ret->y -= params->grid_length;
            depend_extent_ret->h += params->grid_length;
            break;
    }

    printf("Extent initialized for node %d!\ndepend_extent = (%fx%f + %f,%f)\nrecv_extent = (%fx%f + %f,%f)\n",
            rank,
            depend_extent_ret->w, depend_extent_ret->h, depend_extent_ret->x, depend_extent_ret->y,
            recv_extent_ret->w, recv_extent_ret->h, recv_extent_ret->x, recv_extent_ret->y);
}
