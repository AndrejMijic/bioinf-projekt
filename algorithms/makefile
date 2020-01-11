
main.exe : kmer.o hirschberg.o main.cpp 
	g++ -o $@ $^ -O3

kmer.o : kmer.cpp 
	g++ -o $@ $^ -c -O3

hirschberg.o : hirschberg.cpp 
	g++ -o $@ $^ -c -O3


.PHONY: kmer

kmer: kmer.cpp kmer.h
	g++ -o kmer.exe  $^ -DKMER -O3
