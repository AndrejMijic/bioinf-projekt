#ifndef MAIN_H_
#define MAIN_H_

/*
	Contains mutation votes from sequencing result for an index in the genome.
*/
struct occurrence_t {
    int A = 0;
    int C = 0;
    int G = 0;
    int T = 0;
    int del = 0;
    std::vector<struct occurrence_t> insert;
};

#endif // MAIN_H_