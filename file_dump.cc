#include "Particle.h"

#include <fcntl.h>
#include <unistd.h>
#include <string>
#include <stdio.h>
#include <assert.h>

static int fd = -1;

void
dump_init(int rank)
{
    std::string file_name = "nbody_dump_node_" + std::to_string(rank) + ".part.bin";
    fd = open(file_name.c_str(), O_WRONLY | O_CREAT, 0777);
}

void
dump(Particle* part_arr, int n, Simulator_Params* params)
{
    assert(fd != -1);
    double nn = static_cast<double>(n);
    write(fd, &nn, sizeof(double));

    for (int i = 0; i < params->n_particles*4; i++) {
        if (!part_arr[i].invalid) {
            n--;
            write(fd, &part_arr[i].position, sizeof(Vector));
        }
    }
}

void
dump_clean(void)
{
    close(fd);
}
