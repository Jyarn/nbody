EXE=nbody
OBJ=main.o
CC=g++
CXXFLAGS=-g -fopenmp
LDFLAGS=-lraylib

run: $(EXE)
	./$(EXE)

$(EXE): $(OBJ)
	$(CC) $(CXXFLAGS) $(LDFLAGS) $(OBJ) -o $(EXE)

%.o: %.cc
	$(CC) $(CXXFLAGS) -c $? -o $@

clean:
	rm *.o test
