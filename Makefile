EXE=nbody
OBJ=main.o Particle.o equations.o io.o
CC=mpic++
CXXFLAGS=-g
LDFLAGS=-lraylib

build: $(EXE)

debug:
	mpirun -np 4 alacritty -e gdb ./$(EXE)

run: $(EXE)
	mpirun -n 4 ./$(EXE)

$(EXE): $(OBJ)
	$(CC) $(CXXFLAGS) $(LDFLAGS) $(OBJ) -o $(EXE)

%.o: %.cc
	$(CC) $(CXXFLAGS) -c $? -o $@

clean:
	rm *.o test vgcore*
