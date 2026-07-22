MPICXX=mpicxx
MPICXXFLAGS= -cxx=icpc -std=c++11 -O3 -fpermissive -fPIC #-Wall 
INCLUDES= -I../eigen3.37 -I../galib247/include
LDFLAGS= -L../galib247/lib
LIBS= -lga #-mkl

all: erepopt

erepopt: erepopt.o erepobj.o input.o output.o molecule.o skf.o tools.o ga.o genes.o erepobj.hpp molecule.hpp element.hpp tools.hpp ga.hpp ddh.hpp genes.hpp
	$(MPICXX) $(MPICXXFLAGS) $(INCLUDES) -o erepopt erepopt.o erepobj.o input.o output.o molecule.o skf.o tools.o ga.o genes.o $(LDFLAGS) $(LIBS) 

erepopt.o: erepopt.cpp erepobj.hpp molecule.hpp element.hpp tools.hpp ga.hpp ddh.hpp
	$(MPICXX) $(MPICXXFLAGS) $(INCLUDES) -c erepopt.cpp  

erepobj.o: erepobj.cpp erepobj.hpp molecule.hpp element.hpp
	$(MPICXX) $(MPICXXFLAGS) $(INCLUDES) -c erepobj.cpp  

input.o: input.cpp erepobj.hpp molecule.hpp element.hpp ddh.hpp
	$(MPICXX) $(MPICXXFLAGS) $(INCLUDES) -c input.cpp  

output.o: output.cpp erepobj.hpp molecule.hpp element.hpp ddh.hpp
	$(MPICXX) $(MPICXXFLAGS) $(INCLUDES) -c output.cpp  

molecule.o: molecule.cpp erepobj.hpp molecule.hpp element.hpp tools.hpp
	$(MPICXX) $(MPICXXFLAGS) $(INCLUDES) -c molecule.cpp  

skf.o: skf.cpp erepobj.hpp  tools.hpp
	$(MPICXX) $(MPICXXFLAGS) $(INCLUDES) -c skf.cpp  

tools.o:  tools.cpp tools.hpp 
	$(MPICXX) $(MPICXXFLAGS) $(INCLUDES) -c tools.cpp 

ga.o:  ga.cpp ga.hpp erepobj.hpp genes.hpp
	$(MPICXX) $(MPICXXFLAGS) $(INCLUDES) -c ga.cpp 

genes.o: genes.cpp genes.hpp erepobj.hpp ddh.hpp
	$(MPICXX) $(MPICXXFLAGS) $(INCLUDES) -c genes.cpp

clean:
	rm -rf erepopt *.o *.log *.d

