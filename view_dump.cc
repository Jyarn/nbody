#include <unistd.h>
#include <raylib.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <assert.h>
#include <stdio.h>

struct file {
    int size;
    double* data;
    int index;
};

void
load_file(const char* fname, struct file *file_ret)
{
    int fd = open(fname, O_RDONLY);
    assert(fd != -1);

    struct stat s;
    fstat(fd, &s);

    file_ret->size = s.st_size;
    file_ret->index = 0;

    file_ret->data = new double[file_ret->size / 4];
    read(fd, file_ret->data, file_ret->size);
}

int
draw_circles(struct file* file_ret, Color colour)
{
    if (file_ret->size == file_ret->index)
        return 0;

    int index = file_ret->index;
    int n = index + static_cast<int>(file_ret->data[index]) * 2 + 1;
    index++;

    for (; index < n; index += 2) {
        int x = static_cast<double>(file_ret->data[index]*4);
        int y = static_cast<double>(file_ret->data[index+1]*4);
        //DrawCircle(x, y, 10, colour);
    }

    file_ret->index = index;
    assert(file_ret->size >= index);

    return 1;
}

int
main(void)
{
    const char* node_1_name = "nbody_dump_node_0.part.bin";
    const char* node_2_name = "nbody_dump_node_1.part.bin";
    const char* node_3_name = "nbody_dump_node_2.part.bin";
    const char* node_4_name = "nbody_dump_node_3.part.bin";

    struct file dumps[4];
    load_file(node_1_name, &dumps[0]);
    load_file(node_2_name, &dumps[1]);
    load_file(node_3_name, &dumps[2]);
    load_file(node_4_name, &dumps[3]);

    Color colours[4] = { DARKBLUE, DARKBROWN, DARKGREEN, DARKPURPLE };

    //InitWindow(1920, 1080, "nbody");
    //SetTargetFPS(144);
    int run = 1;

    //while (!WindowShouldClose() && run) {
    while(run) {
        //BeginDrawing();
        //ClearBackground(RAYWHITE);
        for (int i = 0; i < 4; i++)
            run = draw_circles(&dumps[i], colours[i]);
        //EndDrawing();
    }

    //CloseWindow();
    return 0;
}
