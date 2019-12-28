#include "kmer.h"
#include "hirschberg.h"
#include <thread>
#include <mutex>

#define PRINT

typedef struct occurrence
{
    int A = 0;
    int C = 0;
    int G = 0;
    int T = 0;
} occurrence_t;

std::mutex file_mtx;
std::mutex occurrences_mtx;

void sequence_to_reference_map(std::ifstream &sequences, std::string &reference,
                               std::map<std::string, int> &minimizer_map,
                               std::vector<int> &minimizers,
                               std::vector<int> &minimizers_locations,
                               std::vector<occurrence_t> &nucliobase_occurrence,
                               int kmer_size, int minimizer_size, int x)
{
    printf("ole\n");
    std::string sequence;
    std::string reverse;
    std::vector<int> sequence_minimizers;
    std::vector<int> sequence_minimizers_reverse;

    printf("%d\n",x);

    while (1)
    {
    
        // thread saftey
        file_mtx.lock();
        int ret = read_sequence(sequences, sequence);
        //end thread saftey
        //printf("%d\n",x);
        file_mtx.unlock();
        
    
        if (ret == 2)
            break; // final sequence was read

        if (ret == 1)
            break; // another thread closed the file
        
        create_complement(reverse, sequence);
        
        minimize_sequence(sequence, kmer_size, minimizer_size, minimizer_map, sequence_minimizers);
        minimize_sequence(reverse, kmer_size, minimizer_size, minimizer_map, sequence_minimizers_reverse);

        subsequence_info_t sequence_rez = subsequence_size<int>(&minimizers[0], minimizers.size(),
                                                    &sequence_minimizers[0], sequence_minimizers.size());

        subsequence_info_t reverse_rez = subsequence_size<int>(&minimizers[0], minimizers.size(),
                                                        &sequence_minimizers_reverse[0],
                                                        sequence_minimizers_reverse.size());

#ifdef PRINT
        printf("%ld %ld %ld %d %d \n",  sequence_rez.length, sequence_rez.first_index,
                                        sequence_rez.last_index, minimizers_locations[sequence_rez.first_index],
                                        minimizers_locations[sequence_rez.last_index]);

        printf("%ld %ld %ld %d %d \n\n",  reverse_rez.length, reverse_rez.first_index,
                                        reverse_rez.last_index, minimizers_locations[reverse_rez.first_index],
                                        minimizers_locations[reverse_rez.last_index]);
#endif
        std::string *better_sequence;
        subsequence_info_t better_rez;

        if (sequence_rez.length >= reverse_rez.length)
        {
            better_sequence = &sequence;
            better_rez = sequence_rez;
        }
        else
        {
            better_sequence = &reverse;
            better_rez = reverse_rez;
        }
        
        std::pair<std::string, std::string> r = needlemanWunsch(reference.substr(minimizers_locations[better_rez.first_index],
                                                                minimizers_locations[better_rez.last_index]),
                                                                *better_sequence, -1, -2, -2, 2);

        // thread saftey
        occurrences_mtx.lock();
        for (int i = 0; i <= better_rez.last_index - better_rez.first_index; i++)
        {
            switch (r.second[i])
            {
            case 'A':
                nucliobase_occurrence[better_rez.first_index + i].A++;
                break;
            case 'C':
                nucliobase_occurrence[better_rez.first_index + i].C++;
                break;
            case 'G':
                nucliobase_occurrence[better_rez.first_index + i].G++;
                break;
            case 'T':
                nucliobase_occurrence[better_rez.first_index + i].T++;
                break;
            default:
                break;
            }
        }
        //end thread saftey
        occurrences_mtx.unlock();
    }
}

int main(int argc, char const *argv[])
{

    int kmer_size = 15, minimizer_size = 5;
    std::string reference;
    std::string sequence;
    char file_line[256];
    std::ifstream file;
    std::ifstream sequences;

    std::map<std::string, int> minimizer_map;
    std::map<int,std::vector<minimizer_info_t>> index_map;
    std::vector<int> minimizers;
    std::vector<int> minimizers_locations;
    std::vector<int> sequence_minimizers;
    std::vector<int> sequence_minimizers_reverse;
    std::string reverse;
    std::string reference_file("data/lambda.fasta");
    std::string sequence_file("data/lambda_simulated_reads.fasta");

    int ret = read_reference(reference_file, reference);
    if (ret == -1)
        return 0;
    
    std::vector<occurrence_t> nucliobase_occurrence(reference.size());

    minimize_reference(reference, kmer_size, minimizer_size, minimizer_map, index_map, minimizers, minimizers_locations);


    sequences.open(sequence_file, std::ifstream::in);


    std::vector<std::thread> threads;

    for (int i = 0; i < 8; i++)
    {
        threads.push_back(std::thread(sequence_to_reference_map, std::ref(sequences),
                                      std::ref(reference), std::ref(minimizer_map), std::ref(minimizers),
                                      std::ref(minimizers_locations),std::ref(nucliobase_occurrence),
                                      kmer_size, minimizer_size, i));
    }

    for ( auto &thread : threads)
    {
        thread.join();
    }



}
