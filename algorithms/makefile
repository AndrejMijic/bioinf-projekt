
main.exe : kmer.o hirschberg.o main.cpp
	g++ -o $@ $^

kmer.o : kmer.cpp 
	g++ -o $@ $^ -c

hirschberg.o : hirschberg.cpp 
	g++ -o $@ $^ -c


.PHONY: kmer

kmer: kmer.cpp kmer.h
	g++ -o kmer.exe  $^ -DKMER
