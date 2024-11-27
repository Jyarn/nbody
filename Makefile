EXE=nbody_
OBJ=main.o Particle.o io.o
CC=mpic++
CXXFLAGS=-O3 -Wall -Wextra -Wpedantic -fopenmp
LDFLAGS=

build: $(EXE)

val: $(EXE)
	mpirun -n 4 valgrind --track-origins=yes --suppressions=/usr/share/openmpi/openmpi-valgrind.supp ./$(EXE)

view_dump: view_dump.cc
	g++ -O3 -lraylib -Wall -Wextra -Wpedantic view_dump.cc -o view_dump
	./view_dump

debug:
	mpirun -np 4 alacritty -e gdb ./$(EXE)

run: $(EXE)
	mpirun -n 4 ./$(EXE)

$(EXE): $(OBJ)
	$(CC) $(CXXFLAGS) $(LDFLAGS) $(OBJ) -o $(EXE)

%.o: %.cc
	$(CC) $(CXXFLAGS) -c $? -o $@

clean:
	rm *.o test vgcore* $(EXE)
