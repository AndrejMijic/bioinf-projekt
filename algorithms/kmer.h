#ifndef _KMER_H
#define _KMER_H

#include <stdio.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <string.h>
#include <map>
#include <vector>
#include <limits>
#include "longest_common_subsequence.h"


typedef struct minimizer_info
{
    int index;
    int num_of_frames;
} minimizer_info_t;

int compare(const char *a, const char *b, const int length);
int find_minimizer(const char *kmer_array, int kmer_size, std::string &minimizer, int minimizer_size);
void minimize_reference(std::string seq, int kmer_size, int minimizer_size,
                        std::map<std::string, int> &minimizer_map, std::map<int,std::vector<minimizer_info_t>> &index_map,
                        std::vector<int> &minimizers,  std::vector<int> &minimizers_locations);
void minimize_sequence(std::string seq, int kmer_size, int minimizer_size,
                       std::map<std::string, int> &minimizer_map,
                       std::vector<int> &minimizers ,std::vector<int> &minimizers_locations);
void create_complement(std::string &complement,std::string original);
int read_reference(std::string& file_path, std::string& reference);
int read_sequence(std::ifstream& sequences_file,std::string &sequence);

#endif //_KMER_H