CC=g++
#CFLAGS=-I.
#CFLAGS=-c -Wall -std=c++0x
CFLAGS=-c -Wall 
LFLAGS=
OPTIMIZE = -O3
#DEPS = fsi.h fluid.h solid.h node.h region.h
#OBJ = main.o fsi.o fluid.o solid.o node.o region.o

DEPS = lb.h boundary.h units.h cell.h ibm.h solid.h chain.h random_mars.h
OBJ = main.o lb.o boundary.o units.o cell.o ibm.o solid.o chain.o random_mars.o
 
%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) $(OPTIMIZE)

fsi: $(OBJ)
	$(CC) -o $@ $^ $(LFLAGS) $(OPTIMIZE)

clean:
	rm -f fsi vel.dat *.o *~
remove:
	rm fluidRst.txt fluidForce.txt cellRst.txt cellForce.txt fgeom.txt cellVelocity.txt chainRst.txt chainForce.txt chainVelocity.txt Log.txt
