CC=g++ -std=c++11
INCLUDE = lattice.hpp observables.hpp iparameters.hpp hamiltonian.hpp runningaverage.hpp
default: Y2-12-7
clean:
	rm -f *.o
# .o files depend on .cpp file of same name + includes
%.o:: %.cpp $(INCLUDE)
#	$(CC) -c -g $< -o $@
	$(CC) -c -O3 $< -o $@
Y2-12-7: Y2-12-7.o lattice.o observables.o iparameters.o hamiltonian.o runningaverage.o
	g++ -std=c++11 $^ -o $@ 
