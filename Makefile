CC=g++ -std=c++11

INCLUDE = lattice.hpp

default: Y2-12-7

clean:
	rm -f *.o

# .o files depend on .cpp file of same name + includes
%.o:: %.cpp $(INCLUDE)
	$(CC) -c $< -o $@

Y2-12-7: Y2-12-7.o lattice.o
	g++ -std=c++11 $^ -o $@
