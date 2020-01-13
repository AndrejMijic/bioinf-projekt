mutation_checker : kmer.o hirschberg.o src/main.cpp output.o
	g++ -std=c++17 -O3 -o $@ $^ -lpthread

kmer.o : src/kmer.cpp 
	g++ -std=c++17 -O3 -o $@ $^ -c -lpthread

hirschberg.o : src/hirschberg.cpp 
	g++ -std=c++17 -O3 -o $@ $^ -c -lpthread 

output.o : src/output.cpp
	g++ -std=c++17 -O3 -o $@ $^ -c

clean:
	rm *.o