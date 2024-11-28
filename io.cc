#include "io.h"

#include <mpi.h>
#include <stdio.h>

static int iter_no = 0;
static Extent depend_extent_arr[4];
static Extent partition_extent_arr[4];

void
send_particle(int index, Particle part, int recv_rank)
{
    int err = MPI_Send(&index, 1, MPI_INT, recv_rank, iter_no, MPI_COMM_WORLD);
    assert(err == MPI_SUCCESS);

    double dbl_buf[SYNC_BUFFER_SIZE];
    dbl_buf[SYNC_POS_X] = part.position.x;
    dbl_buf[SYNC_POS_Y] = part.position.y;
    dbl_buf[SYNC_VEL_X] = part.velocity.x;
    dbl_buf[SYNC_VEL_Y] = part.velocity.y;
    dbl_buf[SYNC_ACC_X] = part.acceleration.x;
    dbl_buf[SYNC_ACC_Y] = part.acceleration.y;
    dbl_buf[SYNC_MASS] = part.mass;

    err = MPI_Send(dbl_buf, SYNC_BUFFER_SIZE, MPI_DOUBLE, recv_rank, iter_no, MPI_COMM_WORLD);
    assert(err == MPI_SUCCESS);
}

void
sync_terminate_connection(int recv_rank)
{
    int term = SYNC_TERM;
    MPI_Ssend(&term, 1, MPI_INT, recv_rank, iter_no, MPI_COMM_WORLD);
}

int
sync_send(Particle* part_arr, int recv_rank, Simulator_Params* params)
{
    Particle part_i;
    Vector part_pos;
    int n_particles_lost = 0;

    Extent partition_extent, depend_extent;
    partition_extent = partition_extent_arr[recv_rank];
    depend_extent = depend_extent_arr[recv_rank];

    for (int i = 0; i < params->n_particles*4; i++) {
        part_i = part_arr[i];
        part_pos = part_i.position;

        // skip all valid and dependency particles
        if (part_i.depend) {
            // assert to make sure we aren't skipping particles just because they have certain flags enabled
//            assert((part_i.depend && (IN_EXTENT_X(part_pos.x, depend_extent) && IN_EXTENT_Y(part_pos.y, depend_extent))) || !part_i.depend);
            continue;
        }

        if (IN_EXTENT_X(part_pos.x, partition_extent) && IN_EXTENT_Y(part_pos.y, partition_extent))
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

        } else if (!part_i.invalid && IN_EXTENT_X(part_pos.x, depend_extent) && IN_EXTENT_Y(part_pos.y, depend_extent)) {
            // potentially unneeded since a particle can be out of the current nodes partition (ie. invalid) but also a dependency particle for other nodes
//            assert(!part_i.invalid);
            assert(IN_EXTENT_X(part_pos.x, partition_extent_arr[params->self_rank]));
            assert(IN_EXTENT_Y(part_pos.y, partition_extent_arr[params->self_rank]));
            send_particle(i, part_i, recv_rank);

        } else {
            //assert((part_pos.x == -1 && part_pos.y == -1) || (IN_EXTENT_X(part_pos.x, partition_extent_arr[params->self_rank]) || IN_EXTENT_Y(part_pos.y, partition_extent_arr[params->self_rank])));
        }

    }

    sync_terminate_connection(recv_rank);
    return n_particles_lost;
}

bool
sync_receive_particle(int sender_rank, int* index_ret, double* dbl_ret)
{
    MPI_Status status;
    int err = MPI_Recv(index_ret, 1, MPI_INT, sender_rank, iter_no, MPI_COMM_WORLD, &status);

    if (*index_ret == SYNC_TERM)
        return false; // return false to terminate loop

    assert(err == MPI_SUCCESS);

    err = MPI_Recv(dbl_ret, SYNC_BUFFER_SIZE, MPI_DOUBLE, sender_rank, iter_no, MPI_COMM_WORLD, &status);
    assert(err == MPI_SUCCESS);

    return true;
}

int
sync_receive(Particle* part_arr, Particles* particle_hash_map, int sender_rank,
        Simulator_Params* params)
{
    int part_index = -1;
    double dbl_buf[SYNC_BUFFER_SIZE] = { 0, 0, 0, 0, 0, 0, 0 };
    int n_particles_gained = 0;

    Extent partition_extent, depend_extent;
    partition_extent = partition_extent_arr[params->self_rank];
    depend_extent = depend_extent_arr[params->self_rank];

    while (true) {
        if (!sync_receive_particle(sender_rank, &part_index, dbl_buf))
            break;

        if (IN_EXTENT_X(dbl_buf[SYNC_POS_X], partition_extent) && IN_EXTENT_Y(dbl_buf[SYNC_POS_Y], partition_extent))
        {
            assert(0 <= part_index && part_index < params->n_particles*4);

            n_particles_gained++;

            part_arr[part_index].depend = false;
            part_arr[part_index].velocity.x = dbl_buf[SYNC_VEL_X];
            part_arr[part_index].velocity.y = dbl_buf[SYNC_VEL_Y];
            part_arr[part_index].acceleration.x = dbl_buf[SYNC_ACC_X];
            part_arr[part_index].acceleration.y = dbl_buf[SYNC_ACC_Y];
            part_arr[part_index].mass = dbl_buf[SYNC_MASS];

            if (part_arr[part_index].invalid) {
                part_arr[part_index].position.x = dbl_buf[SYNC_POS_X];
                part_arr[part_index].position.y = dbl_buf[SYNC_POS_Y];
                insert_particle(particle_hash_map, &part_arr[part_index], params, depend_extent);

            } else
                move_particle(particle_hash_map, &part_arr[part_index], params, depend_extent, dbl_buf[SYNC_POS_X], dbl_buf[SYNC_POS_Y]);

            assert(part_arr[part_index].depend == false && !part_arr[part_index].invalid);

        // if particle is in the dependency extent, sync without deleting the particle
        } else if (IN_EXTENT_X(dbl_buf[SYNC_POS_X], depend_extent) && IN_EXTENT_Y(dbl_buf[SYNC_POS_Y], depend_extent))
        {
            int part_start = sender_rank * params->n_particles;
            int part_len = params->n_particles;

            // assert particle is a valid index in part_arr and part_arr[part_index]
            // this shouldn't hold
            //assert(part_start <= part_index);
            //assert(part_index < part_start + part_len);

            // assert(!part.invalid => part.depend);
            assert((!part_arr[part_index].invalid && part_arr[part_index].depend) || part_arr[part_index].invalid);

            // mark particle as a dependency and not to calculate acceleration for it during one_step
            part_arr[part_index].depend = true;
            part_arr[part_index].velocity.x = dbl_buf[SYNC_VEL_X];
            part_arr[part_index].velocity.y = dbl_buf[SYNC_VEL_Y];
            part_arr[part_index].acceleration.x = dbl_buf[SYNC_ACC_X];
            part_arr[part_index].acceleration.y = dbl_buf[SYNC_ACC_Y];
            part_arr[part_index].mass = dbl_buf[SYNC_MASS];

            if (part_arr[part_index].invalid) {
                part_arr[part_index].position.x = dbl_buf[SYNC_POS_X];
                part_arr[part_index].position.y = dbl_buf[SYNC_POS_Y];
                insert_particle(particle_hash_map, &part_arr[part_index], params, depend_extent);

            } else
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
sync_all(Particle* part_arr, Particles* hash_map, Simulator_Params* params)
{
    int delta_particles = 0;

    switch (params->self_rank) {
        case 0:
            // sync(0, 1)
            delta_particles -= sync_send(part_arr, 1, params);
            delta_particles += sync_receive(part_arr, hash_map, 1, params);

            // sync(0, 2)
            delta_particles -= sync_send(part_arr, 2, params);
            delta_particles += sync_receive(part_arr, hash_map, 2, params);

            // sync(0, 3)
            delta_particles -= sync_send(part_arr, 3, params);
            delta_particles += sync_receive(part_arr, hash_map, 3, params);
            break;

        case 1:
            // sync(1, 0)
            delta_particles += sync_receive(part_arr, hash_map, 0, params);
            delta_particles -= sync_send(part_arr, 0, params);

            // sync(1, 3)
            delta_particles += sync_receive(part_arr, hash_map, 3, params);
            delta_particles -= sync_send(part_arr, 3, params);

            // sync(1, 2)
            delta_particles += sync_receive(part_arr, hash_map, 2, params);
            delta_particles -= sync_send(part_arr, 2, params);
            break;

        case 2:
            // sync(2, 3)
            delta_particles -= sync_send(part_arr, 3, params);
            delta_particles += sync_receive(part_arr, hash_map, 3, params);

            // sync(2, 0)
            delta_particles += sync_receive(part_arr, hash_map, 0, params);
            delta_particles -= sync_send(part_arr, 0, params);

            // sync(2, 1)
            delta_particles -= sync_send(part_arr, 1, params);
            delta_particles += sync_receive(part_arr, hash_map, 1, params);
            break;

        case 3:
            // sync(3, 2)
            delta_particles += sync_receive(part_arr, hash_map, 2, params);
            delta_particles -= sync_send(part_arr, 2, params);

            // sync(3, 1)
            delta_particles -= sync_send(part_arr, 1, params);
            delta_particles += sync_receive(part_arr, hash_map, 1, params);

            // sync(3, 0)
            delta_particles += sync_receive(part_arr, hash_map, 0, params);
            delta_particles -= sync_send(part_arr, 0, params);
            break;
    }

    iter_no++;
    return delta_particles;
}

void
sync_init(Simulator_Params* params, Extent* depend_extent_ret,
        Extent* partition_extent_ret)
{
    assert(params->self_rank < 4);

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

    *depend_extent_ret = depend_extent_arr[params->self_rank];
    *partition_extent_ret = partition_extent_arr[params->self_rank];
}
