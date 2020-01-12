
main.exe : kmer.o hirschberg.o main.cpp output.o
	g++ -std=c++17 -O3 -lpthread -pg -o $@ $^ -O3 -lpthread -pg

kmer.o : kmer.cpp 
	g++ -std=c++17 -O3 -lpthread -pg -o $@ $^ -c -O3 -lpthread -pg

hirschberg.o : hirschberg.cpp 
	g++ -std=c++17 -O3 -lpthread -pg -o $@ $^ -c 


.PHONY: kmer

kmer: kmer.cpp kmer.h
	g++ -std=c++17 -o kmer.exe -DKMER -O3 -lpthread -pg $^ 

output.o : output.cpp
	g++ -std=c++17 -O3 -pg -o $@ $^ -c
