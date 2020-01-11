#include <stdio.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <string.h>
#include <map>
#include <vector>
#include "longest_common_subsequence.h"
#include "longest_common_substring.h"

typedef struct minimizer_info
{
    int index;
    int num_of_frames;
} minimizer_info_t;

int compare(const char *a, const char *b, const int length)
{
    int i = 0;
    for ( i = 0; i < length; i++)
    {
        if (a[i] < b[i])
            return -1;
        else if (a[i] > b[i])
            return 1;
    }

    return 0;
}

int find_minimizer(const char *kmer_array, int kmer_size, std::string &minimizer, int minimizer_size)
{
    int index;
    minimizer.clear();
    for (int i = 0; i < minimizer_size; i++)
    {
        minimizer.push_back('Z');
    }
    for (int i = 0; i < kmer_size - minimizer_size + 1; i++)
    {
        if (compare(kmer_array+i, minimizer.c_str(), minimizer_size) < 0)
        {
            minimizer.clear();
            index = i;
            minimizer.append(kmer_array+i, minimizer_size);
        }
    }

    return index;
}


void minimize_reference(std::string seq, int kmer_size, int minimizer_size,
                       std::map<std::string, int> &minimizer_map, std::map<int,std::vector<minimizer_info_t>> &index_map,
                       std::vector<int> &minimizers )
{
    int unique_id = 0;
    int minimizer_index;
    std::string minimizer;
    std::string prev_minimizer;
    
    const char *sequence = seq.c_str();
    for (int i = 0; i < minimizer_size; i++)
    {
        minimizer.push_back('Z');
        prev_minimizer.push_back('Z');
    }

    for (int i = 0; i < seq.size() - kmer_size + 1; i++)
    {
        std::string x = seq.substr(i + kmer_size - minimizer_size, minimizer_size);
        if (i == 0 || minimizer_index < i)
            minimizer_index = i + find_minimizer(&sequence[i], kmer_size, minimizer, minimizer_size);
        else if (x.compare(minimizer) < 0)
        {
            minimizer_index = i + kmer_size - minimizer_size;
            minimizer.clear();
            minimizer.assign(seq, i + kmer_size - minimizer_size, minimizer_size);
        }

        //std::cout << f << " " << minimizer.c_str() << "  " << prev_minimizer.c_str() << std::endl;
        
        if (minimizer.compare(prev_minimizer) != 0)
        {
            std::map<std::string, int>::iterator it =  minimizer_map.find(minimizer);
            if ( it ==  minimizer_map.end())
            {
                minimizer_map[minimizer] = unique_id;
                std::vector<minimizer_info_t> vec;
                minimizer_info_t info = {i, 1};
                vec.push_back(info);
                index_map[unique_id] = vec;
                minimizers.push_back(unique_id);
                unique_id++;
            }
            else
            {
                int id = it->second;
                minimizer_info_t info = {i, 1};
                index_map[id].push_back(info);
                minimizers.push_back(id);
            }
        }
        
        else
        {
            std::map<std::string, int>::iterator it =  minimizer_map.find(minimizer);
            int id = it->second;
            std::vector<minimizer_info_t> vec = index_map[id];
            vec.back().num_of_frames++;
            minimizers.push_back(id);
        }
        prev_minimizer.clear();
        prev_minimizer.assign(minimizer);
        
    }
    
}

void minimize_sequence(std::string seq, int kmer_size, int minimizer_size,
                       std::map<std::string, int> &minimizer_map, std::map<int,std::vector<minimizer_info_t>> &index_map,
                       std::vector<int> &minimizers )
{
    int minimizer_index;
    int unknown_minimizer_index = -1;
    std::string minimizer;
    std::string prev_minimizer;

    const char *sequence = seq.c_str();
    for (int i = 0; i < minimizer_size; i++)
    {
        minimizer.push_back('Z');
        prev_minimizer.push_back('Z');
    }

    for (int i = 0; i < seq.size() - kmer_size + 1; i++)
    {
        std::string x = seq.substr(i + kmer_size - minimizer_size, minimizer_size);
        if (i == 0 || minimizer_index < i)
            minimizer_index = i + find_minimizer(&sequence[i], kmer_size, minimizer, minimizer_size);
        else if (x.compare(minimizer) < 0)
        {
            minimizer_index = i + kmer_size - minimizer_size;
            minimizer.clear();
            minimizer.assign(seq, i + kmer_size - minimizer_size, minimizer_size);
        }

        std::map<std::string, int>::iterator it = minimizer_map.find(minimizer);
        if (it == minimizer_map.end())
        {
            minimizers.push_back(unknown_minimizer_index);
            unknown_minimizer_index--;
        }
        else
        {
            int id = it->second;
            minimizers.push_back(id);
        }

        //std::cout <<minimizer <<"  "<<minimizers.back()<<std::endl;

        prev_minimizer.clear();
        prev_minimizer.assign(minimizer);
    }
}

void create_complement(std::string &complement,std::string original)
{
    std::reverse(original.begin(), original.end());

    for (int i = 0; i < original.size(); i++)
    {
        if (original[i] == 'A')
            complement.push_back('T');
        else if (original[i] == 'T')
            complement.push_back('A');
        else if (original[i] == 'G')
            complement.push_back('C');
        else if (original[i] == 'C')
            complement.push_back('G');
    }
}



int main()
{
    int kmer_size = 10, minimizer_size = 4;
    std::string reference = "ATTACGCGATACGTAGCATGCGTAGCGTATTTACGTACAACGTGAACGTGTCAAAGTCCTTCACA";
    std::string sequence = "GCGTAACGATTTACGTACA"; //original GCGTAGCGTATTTACGTACA
    char file_line[256];
    std::ifstream file;
    std::ifstream sequences;

    std::map<std::string, int> minimizer_map;
    std::map<int,std::vector<minimizer_info_t>> index_map;
    std::vector<int> minimizers;
    std::vector<int> sequence_minimizers;

    minimize_reference(reference, kmer_size, minimizer_size, minimizer_map, index_map, minimizers);
    minimize_sequence(sequence, kmer_size, minimizer_size, minimizer_map, index_map, sequence_minimizers);


    int *a = &minimizers[0];
    int *b = &sequence_minimizers[0];

    for(auto el : minimizers)
    {
        printf("%2d ",el);
    }
    printf("\n");

    for(auto el : sequence_minimizers)
    {
        printf("%2d ",el);
    }
    printf("\n");


    subsequence_info_t rez = subsequence_size<int>(a, minimizers.size(), b, sequence_minimizers.size());
    printf("Subsequence\n Reference first index: %ld, last_index: %ld, length: %ld\n", rez.first_index_ref, rez.last_index_ref, rez.length);

    substring_info_t rez2 = substring_size<int>(a, minimizers.size(), b, sequence_minimizers.size());
    printf("Substring\n Reference first index: %ld, length: %ld\n", rez2.first_index, rez2.length);

    return 0;

    
/*  // reading from file

    file.open("data/lambda.fasta", std::ifstream::in);

    if (!file.is_open()){
        printf("Error opening file.");
        return 1;
    }

    file.getline(file_line, 256);

    std::cout << file_line << "\n";
    reference.clear();
    seqence.clear();
    while (true)
    {
        file.getline(file_line, 256);
        if (strlen(file_line) == 0)
            break;
        kmer.append(file_line);

    }
    minimize_reference(reference, kmer_size, minimizer_size, minimizer_map, index_map, minimizers);

    sequences.open("data/lambda_simulated_reads.fasta", std::ifstream::in);

    if (!sequences.is_open()){
        printf("Error opening file.");
        return 1;
    }

    sequences.getline(file_line, 256);
    
    while(true)
    {
        char c = (char) sequences.get();
        if (c == '\n')
            break;
        sequence.push_back(c);
    }

    minimize_sequence(sequence, kmer_size, minimizer_size, minimizer_map, index_map, sequence_minimizers);

    int *a = &minimizers[0];
    int *b = &sequence_minimizers[0];

    subsequence_info_t rez = subsequence_size<int>(a, minimizers.size(), b, sequence_minimizers.size());

    printf("\n%ld %ld %ld\n",  rez.length, rez.first_index, rez.last_index);


    substring_info_t rez2 = substring_size<int>(a, minimizers.size(), b, sequence_minimizers.size());

    printf("\n%ld %ld %ld\n",  rez2.length, rez2.first_index, rez2.first_index + rez2.length);
*/
}
