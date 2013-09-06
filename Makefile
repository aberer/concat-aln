


CXXFLAGS=-std=c++0x -O3  -Wall -Wextra 

all : concat-fasta

concat-fasta :  concat.o
	$(CXX) $(CXXFLAGS)  -o $@ $< 


clean : 
	rm -f concat-fasta *.o 
