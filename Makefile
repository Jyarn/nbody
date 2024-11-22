EXE=nbody
OBJ=main.o Particle.o equations.o
CC=g++
CXXFLAGS=-g -fopenmp
LDFLAGS=-lraylib

run: $(EXE)
	./$(EXE)

build: $(EXE)

$(EXE): $(OBJ)
	$(CC) $(CXXFLAGS) $(LDFLAGS) $(OBJ) -o $(EXE)

%.o: %.cc
	$(CC) $(CXXFLAGS) -c $? -o $@

clean:
	rm *.o test vgcore*
