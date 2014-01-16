PROG = XXZ
CXX = g++
CXXFLAGS = -Wall -Wextra -O3 -std=c++11 -march=native -I ../Eigen_3.2.0
OBJS = EffectiveHamiltonian.o ESolver.o FreeFunctions.o main.o TheBlock.o XXZ.o
COMMONHS = d.h main.h Hamiltonian.h
light = rm -f *.cpp~ *.h~ Makefile~
deep = rm -f $(PROG) *.o ./Output/*

$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(PROG) $(OBJS)

EffectiveHamiltonian.o: $(COMMONHS) TheBlock.h EffectiveHamiltonian.h FreeFunctions.h ESolver.h

ESolver.o: $(COMMONHS) TheBlock.h ESolver.h

FreeFunctions.o: $(COMMONHS) TheBlock.h EffectiveHamiltonian.h

main.o: $(COMMONHS) TheBlock.h EffectiveHamiltonian.h ESolver.h FreeFunctions.h

TheBlock.o: $(COMMONHS) TheBlock.h EffectiveHamiltonian.h FreeFunctions.h ESolver.h

XXZ.o: $(COMMONHS) 

lightclean:
	$(light)

deepclean:
	$(deep)

clean:
	$(light)
	$(deep)

upload:
	scp *.cpp *.h knot.cnsi.ucsb.edu:~/QTIM
