PROG = XXZ
CXX = g++
CXXFLAGS = -Wall -Wextra -O3 -std=c++11 -march=native -I ../Eigen_3.2.0 $(DEBUG)
LIBS = -llapack
OBJS = EffectiveHamiltonian.o ESolver.o FreeFunctions.o Lanczos.o main.o modifyHamParams.o TheBlock.o XXZ.o
COMMONHS = d.h main.h Hamiltonian.h
light = rm -f *.cpp~ *.h~ Makefile~
git = rm -f $(PROG) ./Output/*
deep = $(git) *.o

$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) $(LIBS) -o $(PROG) $(OBJS)

EffectiveHamiltonian.o: $(COMMONHS) TheBlock.h EffectiveHamiltonian.h FreeFunctions.h ESolver.h

ESolver.o: $(COMMONHS) TheBlock.h ESolver.h

FreeFunctions.o: $(COMMONHS) TheBlock.h EffectiveHamiltonian.h

Lanczos.o: d.h main.h ESolver.h

main.o: $(COMMONHS) TheBlock.h EffectiveHamiltonian.h ESolver.h FreeFunctions.h

TheBlock.o: $(COMMONHS) TheBlock.h EffectiveHamiltonian.h FreeFunctions.h ESolver.h

XXZ.o: $(COMMONHS) 

lightclean:
	$(light)

gitclean:
	$(git)

deepclean:
	$(deep)

clean:
	$(light)
	$(deep)

upload:
	scp *.cpp *.h Makefile knot.cnsi.ucsb.edu:~/XXZ
