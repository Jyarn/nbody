#include "io.h"

#include <assert.h>
#include <mpi.h>

static int iter_no = 0;
static Extent depend_extent_arr[4];
static Extent partition_extent_arr[4];

void
send_particle(int index, Particle part, int recv_rank)
{
    int err = MPI_Send(&index, 1, MPI_INT, recv_rank, iter_no, MPI_COMM_WORLD);
    assert(err == MPI_SUCCESS);

    double dbl_buf[6];
    dbl_buf[SYNC_POS_X] = part.position.x;
    dbl_buf[SYNC_POS_Y] = part.position.y;
    dbl_buf[SYNC_VEL_X] = part.velocity.x;
    dbl_buf[SYNC_VEL_Y] = part.velocity.y;
    dbl_buf[SYNC_ACC_X] = part.acceleration.x;
    dbl_buf[SYNC_ACC_Y] = part.acceleration.y;

    err = MPI_Send(dbl_buf, 6, MPI_DOUBLE, recv_rank, iter_no, MPI_COMM_WORLD);
    assert(err == MPI_SUCCESS);
}

void
sync_terminate_connection(int recv_rank)
{
    int term = SYNC_TERM;
    MPI_Ssend(&term, 1, MPI_INT, recv_rank, iter_no, MPI_COMM_WORLD);
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
            assert((part_i.invalid && !(IN_EXTENT_X(part_pos.x, recv_extent) && IN_EXTENT_Y(part_pos.y, recv_extent))) || !part_i.invalid);
            assert((part_i.depend && (IN_EXTENT_X(part_pos.x, depend_extent) && IN_EXTENT_Y(part_pos.y, depend_extent))) || !part_i.depend);
            continue;
        }

        if (IN_EXTENT_X(part_pos.x, recv_extent) && IN_EXTENT_Y(part_pos.y, recv_extent))
        {
            assert(part_i.invalid);
            assert(!part_i.depend);

            // technically the particle is stil valid (ie. still can send the
            // data in further iteration) thus we should permanently destroy
            // the data so that this isn't possible
            send_particle(i, part_i, recv_rank);
            part_arr[i].position.x = -1;
            part_arr[i].position.y = -1;

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
    int err = MPI_Recv(index_ret, 1, MPI_INT, sender_rank, iter_no, MPI_COMM_WORLD, &status);

    if (*index_ret == SYNC_TERM)
        return false; // return false to terminate loop

    assert(err == MPI_SUCCESS);

    err = MPI_Recv(dbl_ret, 6, MPI_DOUBLE, sender_rank, iter_no, MPI_COMM_WORLD, &status);
    assert(err == MPI_SUCCESS);

    return true;
}

int
sync_receive(Particle* part_arr, int particle_partition_start,
        int particle_partition_len, Particles* particle_hash_map,
        int sender_rank, Extent depend_extent, Extent recv_extent,
        Simulator_Params* params)
{
    int part_index = -1;
    double dbl_buf[6] = { 0, 0, 0, 0, 0, 0 };
    int n_particles_gained = 0;

    while (true) {
        if (!sync_receive_particle(part_arr, sender_rank, &part_index, dbl_buf))
            break;

        if (IN_EXTENT_X(dbl_buf[SYNC_POS_X], recv_extent) && IN_EXTENT_Y(dbl_buf[SYNC_POS_Y], recv_extent))
        {
            assert(0 <= part_index && part_index < params->n_particles*4);
            assert(part_arr[part_index].invalid);

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
            // equivalent to assert(false) but prints an error message
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
    Extent d1, d2, d3, d4;
    Extent r1, r2, r3, r4;
    sync_init(0, params, &d1, &r1);
    sync_init(1, params, &d2, &r2);
    sync_init(2, params, &d3, &r3);
    sync_init(3, params, &d4, &r4);

    int delta_particles = 0;

    switch (rank) {
        case 0:
            // sync(0, 1)
            delta_particles -= sync_send(part_arr, params->n_particles, 1,
                    d2, r2);
            delta_particles += sync_receive(part_arr, particle_partition_start,
                    particle_partiton_len, hash_map, 1, d1,
                    r1, params);

            // sync(0, 2)
            delta_particles -= sync_send(part_arr, params->n_particles, 2,
                    d3, r3);
            delta_particles += sync_receive(part_arr, particle_partition_start,
                    particle_partiton_len, hash_map, 2, d1,
                    r1, params);

            // sync(0, 3)
            delta_particles -= sync_send(part_arr, params->n_particles, 3,
                    d4, r4);
            delta_particles += sync_receive(part_arr, particle_partition_start,
                    particle_partiton_len, hash_map, 3, d1,
                    r1, params);
            break;
        case 1:
            // sync(1, 0)
            delta_particles += sync_receive(part_arr, particle_partition_start,
                    particle_partiton_len, hash_map, 0, d2,
                    r2, params);
            delta_particles -= sync_send(part_arr, params->n_particles, 0,
                    d1, r1);

            // sync(1, 3)
            delta_particles += sync_receive(part_arr, particle_partition_start,
                    particle_partiton_len, hash_map, 3, d2,
                    r2, params);
            delta_particles -= sync_send(part_arr, params->n_particles, 3,
                    d4, r4);

            // sync(1, 2)
            delta_particles += sync_receive(part_arr, particle_partition_start,
                    particle_partiton_len, hash_map, 2, d2,
                    r2, params);
            delta_particles -= sync_send(part_arr, params->n_particles, 2,
                    d3, r3);
            break;

        case 2:
            // sync(2, 3)
            delta_particles -= sync_send(part_arr, params->n_particles, 3,
                    d4, r4);
            delta_particles += sync_receive(part_arr, particle_partition_start,
                    particle_partiton_len, hash_map, 3, d3,
                    r3, params);

            // sync(2, 0)
            delta_particles += sync_receive(part_arr, particle_partition_start,
                    particle_partiton_len, hash_map, 0, d3,
                    r3, params);
            delta_particles -= sync_send(part_arr, params->n_particles, 0,
                    d1, r1);

            // sync(2, 1)
            delta_particles -= sync_send(part_arr, params->n_particles, 1,
                    d2, r2);
            delta_particles += sync_receive(part_arr, particle_partition_start,
                    particle_partiton_len, hash_map, 1, d3,
                    r3, params);
            break;

        case 3:
            // sync(2, 3)
            delta_particles += sync_receive(part_arr, particle_partition_start,
                    particle_partiton_len, hash_map, 2, d4,
                    r4, params);
            delta_particles -= sync_send(part_arr, params->n_particles, 2,
                    d3, r3);

            // sync(3, 1)
            delta_particles -= sync_send(part_arr, params->n_particles, 1,
                    d2, r2);
            delta_particles += sync_receive(part_arr, particle_partition_start,
                    particle_partiton_len, hash_map, 1, d4,
                    r4, params);

            // sync(3, 0)
            delta_particles += sync_receive(part_arr, particle_partition_start,
                    particle_partiton_len, hash_map, 0, d4,
                    r4, params);
            delta_particles -= sync_send(part_arr, params->n_particles, 0,
                    d1, r1);
            break;
    }

    iter_no++;
    return delta_particles;
}

void
sync_init(int rank, Simulator_Params* params, Extent* depend_extent_ret, Extent* partition_extent_ret)
{
    assert(rank < 4);
    partition_extent_arr[0].x = 0.0;
    partition_extent_arr[0].y = 0.0;
    partition_extent_arr[0].w = params->n_cells_x * params->grid_length / 2;
    partition_extent_arr[0].h = params->n_cells_y * params->grid_length / 2;
    depend_extent_arr[0] = partition_extent_arr[0];
    depend_extent_arr[0].w += params->grid_length;
    depend_extent_arr[0].h += params->grid_length;

    partition_extent_arr[1].x = 0.0;
    partition_extent_arr[1].y = params->n_cells_y * params->grid_length / 2;
    partition_extent_arr[1].w = params->n_cells_x * params->grid_length / 2;
    partition_extent_arr[1].h = params->n_cells_y * params->grid_length / 2;
    depend_extent_arr[1] = partition_extent_arr[1];
    depend_extent_arr[1].y -= params->grid_length;
    depend_extent_arr[1].h += params->grid_length;
    depend_extent_arr[1].w += params->grid_length;

    partition_extent_arr[2].x = params->n_cells_x * params->grid_length / 2;
    partition_extent_arr[2].y = 0.0;
    partition_extent_arr[2].w = params->n_cells_x * params->grid_length / 2;
    partition_extent_arr[2].h = params->n_cells_y * params->grid_length / 2;
    depend_extent_arr[2] = partition_extent_arr[2];
    depend_extent_arr[2].x -= params->grid_length;
    depend_extent_arr[2].w += params->grid_length;
    depend_extent_arr[2].h += params->grid_length;

    partition_extent_arr[3].x = params->n_cells_x * params->grid_length / 2;
    partition_extent_arr[3].y = params->n_cells_y * params->grid_length / 2;
    partition_extent_arr[3].w = params->n_cells_x * params->grid_length / 2;
    partition_extent_arr[3].h = params->n_cells_y * params->grid_length / 2;
    depend_extent_arr[3] = partition_extent_arr[3];
    depend_extent_arr[3].x -= params->grid_length;
    depend_extent_arr[3].w += params->grid_length;
    depend_extent_arr[3].y -= params->grid_length;
    depend_extent_arr[3].h += params->grid_length;

    *depend_extent_ret = depend_extent_arr[rank];
    *partition_extent_ret = partition_extent_arr[rank];
}
