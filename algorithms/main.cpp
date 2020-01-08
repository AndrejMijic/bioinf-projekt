#include "kmer.h"
#include "hirschberg.h"
#include <thread>
#include <mutex>

#define PRINT

typedef struct occurrence {
    int A = 0;
    int C = 0;
    int G = 0;
    int T = 0;
} occurrence_t;

std::mutex file_mtx;
std::mutex occurrences_mtx;

//doesn't work but should count what bases occur on what indexes
void update_occurrences(int start_index, int len, std::pair<std::string, std::string> &r,
                       std::vector<occurrence_t> &nucliobase_occurrence) {
    // thread saftey
    occurrences_mtx.lock();
    for (int i = 0; i <= len; i++) {
        switch (r.second[i]) {
            case 'A':
                nucliobase_occurrence[start_index + i].A++;
                break;
            case 'C':
                nucliobase_occurrence[start_index + i].C++;
                break;
            case 'G':
                nucliobase_occurrence[start_index + i].G++;
                break;
            case 'T':
                nucliobase_occurrence[start_index + i].T++;
                break;
            default:
                break;
        }
    }
    //end thread saftey
    occurrences_mtx.unlock();
}

void sequence_to_reference_map(std::ifstream &sequences, std::string &reference,
                               std::map<std::string, int> &minimizer_map,
                               std::vector<int> &minimizers,
                               std::vector<int> &minimizers_locations,
                               std::vector<occurrence_t> &nucliobase_occurrence,
                               int kmer_size, int minimizer_size, int x) {
    std::string sequence;
    std::string reverse;
    std::vector<int> sequence_minimizers;
    std::vector<int> sequence_minimizers_reverse;

    while (1) {
        sequence.clear();
        reverse.clear();
        sequence_minimizers.clear();
        sequence_minimizers_reverse.clear();

        // thread saftey
        file_mtx.lock();
        int ret = read_sequence(sequences, sequence);
        //end thread saftey
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
        if (sequence_rez.length != 0)
        {
            setbuf(stdout, 0);
            printf("%ld %ld %ld %d %d %d\n", sequence_rez.length, sequence_rez.first_index,
                   sequence_rez.last_index, minimizers_locations[sequence_rez.first_index],
                   minimizers_locations[sequence_rez.last_index] + kmer_size,
                   minimizers_locations[sequence_rez.last_index] - minimizers_locations[sequence_rez.first_index] +
                   kmer_size);
        }
        if (reverse_rez.length != 0)
        {
            setbuf(stdout, 0);
            printf("%ld %ld %ld %d %d %d\n\n", reverse_rez.length, reverse_rez.first_index,
                   reverse_rez.last_index, minimizers_locations[reverse_rez.first_index],
                   minimizers_locations[reverse_rez.last_index] + kmer_size,
                   minimizers_locations[reverse_rez.last_index] - minimizers_locations[reverse_rez.first_index] +
                   kmer_size);
        }
#endif
        std::string *better_sequence;
        subsequence_info_t better_rez;

        if (sequence_rez.length >= reverse_rez.length) {

            better_sequence = &sequence;
            better_rez = sequence_rez;
        } else {
            better_sequence = &reverse;
            better_rez = reverse_rez;
        }
        
        std::pair<std::string, std::string> r;

        if (better_rez.length == 0 || sequence.size() * 1.5 <
            minimizers_locations[better_rez.last_index] - minimizers_locations[better_rez.first_index] + kmer_size)
            continue;

        // used because bad_alloc otherwise
        if (minimizers_locations[better_rez.last_index] -
            minimizers_locations[better_rez.first_index] + kmer_size > 10000 && better_sequence->size() > 10000)
        {
            int start_index = 0;
            int len = 10000;
     
            while(true)
            {
                if (start_index >= minimizers_locations[better_rez.last_index] -
                                   minimizers_locations[better_rez.first_index] + kmer_size)
                    break;
                r = needlemanWunsch(
                        reference.substr(minimizers_locations[better_rez.first_index] + start_index, len),
                        better_sequence->substr(start_index, len), -1, -2, -2, 2);

                update_occurrences(start_index, len, r, nucliobase_occurrence);
                start_index = len + start_index;
                if (start_index + len > minimizers_locations[better_rez.last_index] -
                                        minimizers_locations[better_rez.first_index] + kmer_size)
                    len = minimizers_locations[better_rez.last_index] -
                          minimizers_locations[better_rez.first_index] + kmer_size - start_index;
            }
        }
        else
        {
            r = needlemanWunsch(
                    reference.substr(minimizers_locations[better_rez.first_index],
                                     minimizers_locations[better_rez.last_index] -
                                     minimizers_locations[better_rez.first_index] + kmer_size),
                    *better_sequence, -1, -2, -2, 2);

            update_occurrences(minimizers_locations[better_rez.first_index],
                               minimizers_locations[better_rez.last_index] -
                               minimizers_locations[better_rez.first_index] + kmer_size, r, nucliobase_occurrence);
        }


    }
}

int main(int argc, char const *argv[]) {

    int kmer_size = 25, minimizer_size = 10;
    std::string reference;
    std::string sequence;
    char file_line[256];
    std::ifstream file;
    std::ifstream sequences;
    std::string reverse;
    std::vector<int> sequence_minimizers;
    std::vector<int> sequence_minimizers_reverse;

    std::map<std::string, int> minimizer_map;
    std::map<int, std::vector<minimizer_info_t>> index_map;
    std::vector<int> minimizers;
    std::vector<int> minimizers_locations;
    std::string reference_file("data/lambda.fasta");
    std::string sequence_file("data/lambda_simulated_reads.fasta");

    read_reference(reference_file, reference);
    sequences.open(sequence_file, std::ifstream::in);
    minimize_reference(reference, kmer_size, minimizer_size, minimizer_map, index_map, minimizers, minimizers_locations);

    std::vector<occurrence_t> nucliobase_occurrence(reference.size());

    std::vector<std::thread> threads;

    for (int i = 0; i < 8; i++) {
        threads.push_back(std::thread(sequence_to_reference_map, std::ref(sequences),
                                      std::ref(reference), std::ref(minimizer_map), std::ref(minimizers),
                                      std::ref(minimizers_locations), std::ref(nucliobase_occurrence),
                                      kmer_size, minimizer_size, i));
    }

    for (auto &thread : threads) {
        thread.join();
    }

    for (int i = 0; i < reference.size(); i++) {
        occurrence_t occ = nucliobase_occurrence[i];
        printf("%c  %d %d %d %d\n", reference[i], occ.A, occ.C, occ.G, occ.T);
    }


}
