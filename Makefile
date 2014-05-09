CC=g++
CFLAGS=-I.
OPTIMIZE = -O3
#DEPS = fsi.h fluid.h solid.h node.h region.h
#OBJ = main.o fsi.o fluid.o solid.o node.o region.o

DEPS = lb.h boundary.h units.h 
OBJ = main.o lb.o boundary.o units.o
 
%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(OPTIMIZE)

fsi: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(OPTIMIZE)

clean:
	rm -f fsi vel.dat *.o *~
