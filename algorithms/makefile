mutation_checker : kmer.o hirschberg.o main.cpp output.o
	g++ -std=c++17 -O3 -o $@ $^ -lpthread

kmer.o : kmer.cpp 
	g++ -std=c++17 -O3 -o $@ $^ -c -lpthread

hirschberg.o : hirschberg.cpp 
	g++ -std=c++17 -O3 -o $@ $^ -c -lpthread 

output.o : output.cpp
	g++ -std=c++17 -O3 -o $@ $^ -c
