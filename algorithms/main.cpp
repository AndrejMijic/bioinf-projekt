#include "kmer.h"
#include "hirschberg.h"
#include <thread>
#include <mutex>

#define RPRINT

typedef struct occurrence {
    int A = 0;
    int C = 0;
    int G = 0;
    int T = 0;
    int del = 0;
    std::vector<struct occurrence> insert;
} occurrence_t;

std::mutex file_mtx;
std::mutex occurrences_mtx;

//doesn't work but should count what bases occur on what indexes
void update_occurrences(int start_index, int len, std::pair<std::string, std::string> &r,
                       std::vector<occurrence_t> &nucliobase_occurrence)
{

    int index = start_index - 1;
    int insert_index = 0;
    int i = 0;

    if (r.first[i]=='-' )
        while (r.first[i]=='-' ) i++;
    else if (r.second[i]=='-' )
    {
        while (r.second[i]=='-' ) 
        {
            i++;
            index++;
        }
    }

    // thread saftey
    occurrences_mtx.lock();
    for (; i <= r.first.size(); i++) 
    {
        if (r.first[i] == '-')
        {
            if (nucliobase_occurrence[index].insert.size() == insert_index)
            {
                occurrence_t occ;
                nucliobase_occurrence[index].insert.push_back(occ);
            }

            if (r.second[i] == 'A')
            {
                nucliobase_occurrence[index].insert[insert_index].A++;
            }
            else if (r.second[i] == 'C')
            {
                nucliobase_occurrence[index].insert[insert_index].C++;
            }
            else if (r.second[i] == 'G')
            {
                nucliobase_occurrence[index].insert[insert_index].G++;
            }
            else if (r.second[i] == 'T')
            {
                nucliobase_occurrence[index].insert[insert_index].T++;
            }
            insert_index++;

            continue;
        }

        index++;

        insert_index = 0;
        if (r.second[i] == 'A')
        {
            nucliobase_occurrence[index].A++;
        }
        else if (r.second[i] == 'C')
        {
            nucliobase_occurrence[index].C++;
        }
        else if (r.second[i] == 'G')
        {
            nucliobase_occurrence[index].G++;
        }
        else if (r.second[i] == 'T')
        {
            nucliobase_occurrence[index].T++;
        }
        else if (r.second[i] == '-')
        {
            nucliobase_occurrence[index].del++;
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
                               int kmer_size, int minimizer_size, int thread_id) {
    std::string sequence;
    std::string reverse;
    std::vector<int> sequence_minimizers;
    std::vector<int> sequence_minimizers_reverse;
    std::vector<int> minimizers_locations_seq;
    std::vector<int> minimizers_locations_rev;
    std::pair<std::string, std::string> r;
    std::string *better_sequence;
    subsequence_info_t better_rez;
    std::vector<int> *better_minimizer_locations;
    while (1) {
        sequence.clear();
        reverse.clear();
        sequence_minimizers.clear();
        sequence_minimizers_reverse.clear();
        r.first.clear();
        r.second.clear();
        minimizers_locations_seq.clear();
        minimizers_locations_rev.clear();

        // thread saftey
        file_mtx.lock();
        int ret = read_sequence(sequences, sequence);
        //end thread saftey
        file_mtx.unlock();

        if (ret == -2)
            break; // final sequence was read

        if (ret == -1)
            break; // another thread closed the file

        create_complement(reverse, sequence);

        minimize_sequence(sequence, kmer_size, minimizer_size, minimizer_map, sequence_minimizers, minimizers_locations_seq);
        minimize_sequence(reverse, kmer_size, minimizer_size, minimizer_map, sequence_minimizers_reverse, minimizers_locations_rev);

        subsequence_info_t sequence_rez = subsequence_size<int>(&minimizers[0], minimizers.size(),
                                                                &sequence_minimizers[0], sequence_minimizers.size());

        subsequence_info_t reverse_rez = subsequence_size<int>(&minimizers[0], minimizers.size(),
                                                               &sequence_minimizers_reverse[0],
                                                               sequence_minimizers_reverse.size());

#ifdef PRINT
        if (sequence_rez.length != 0)
        {
            setbuf(stdout, 0);
            printf("%d %d %d %d %d %d\n", sequence_rez.length, sequence_rez.first_index,
                   sequence_rez.last_index, minimizers_locations[sequence_rez.first_index],
                   minimizers_locations[sequence_rez.last_index] + kmer_size,
                   minimizers_locations[sequence_rez.last_index] - minimizers_locations[sequence_rez.first_index] +
                   kmer_size);
        }
        if (reverse_rez.length != 0)
        {
            setbuf(stdout, 0);
            printf("%d %d %d %d %d %d\n\n", reverse_rez.length, reverse_rez.first_index,
                   reverse_rez.last_index, minimizers_locations[reverse_rez.first_index],
                   minimizers_locations[reverse_rez.last_index] + kmer_size,
                   minimizers_locations[reverse_rez.last_index] - minimizers_locations[reverse_rez.first_index] +
                   kmer_size);
        }
#endif


        if (sequence_rez.length >= reverse_rez.length) 
        {
            better_sequence = &sequence;
            better_minimizer_locations = &minimizers_locations_seq;
            better_rez = sequence_rez;
        }
        else
        {
            better_sequence = &reverse;
            better_rez = reverse_rez;
            better_minimizer_locations = &minimizers_locations_rev;
        }


        printf("Thread: %d Sequence number: %d Sequence size:%5ld, Reference start: %5d, Reference length: %5d ", thread_id, ret, sequence.size(), minimizers_locations[better_rez.first_index_ref],
        minimizers_locations[better_rez.last_index_ref] - minimizers_locations[better_rez.first_index_ref] + kmer_size);

        if (better_rez.length == 0 || sequence.size() * 2 <
            minimizers_locations[better_rez.last_index_ref] - minimizers_locations[better_rez.first_index_ref] + kmer_size)
        {
            printf("Accepted: NO\n");
            continue;
        }

        printf("Accepted: YES\n");

        r = needlemanWunsch(
                reference.substr(minimizers_locations[better_rez.first_index_ref],
                                    minimizers_locations[better_rez.last_index_ref] -
                                    minimizers_locations[better_rez.first_index_ref] + kmer_size),
                better_sequence->substr(better_minimizer_locations->at(better_rez.first_index_seq),
                better_minimizer_locations->at(better_rez.last_index_seq) -
                better_minimizer_locations->at(better_rez.first_index_seq) + kmer_size), -3, -3, -3, 2);
        //std::cout << r.first << "\n" << r.second << "\n\n";
        update_occurrences(minimizers_locations[better_rez.first_index_ref],
                            minimizers_locations[better_rez.last_index_ref] -
                            minimizers_locations[better_rez.first_index_ref] + kmer_size, r, nucliobase_occurrence);
    }
}

int main(int argc, char const *argv[]) {

    int kmer_size = 26, minimizer_size = 12, num_of_threads = 4;
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
    std::string reference_file;
    std::string sequence_file;
    const char reference_file_c_str[] = "data/lambda.fasta";
    const char sequence_file_c_str[] = "data/lambda_simulated_reads.fasta";
    if (argc == 3)
    {
        kmer_size = atoi(argv[1]);
        minimizer_size = atoi(argv[2]);
    }

    //add getopts if exclusively run on linux
    kmer_size = argc >= 2 ? atoi(argv[1]) : kmer_size;
    minimizer_size = argc >= 3 ? atoi(argv[2]) : minimizer_size;
    num_of_threads = argc >= 4 ? atoi(argv[3]) : num_of_threads;
    reference_file.assign(argc >= 5 ? argv[4] : reference_file_c_str);
    sequence_file.assign(argc >= 6 ? argv[5] : sequence_file_c_str);

    printf("Using kmer size %d, minimizers size %d, %d threads, reference file %s, sequences file %s\n", kmer_size, minimizer_size, num_of_threads, reference_file.c_str(), sequence_file.c_str());

    read_reference(reference_file, reference);
    sequences.open(sequence_file, std::ifstream::in);
    minimize_reference(reference, kmer_size, minimizer_size, minimizer_map, index_map, minimizers, minimizers_locations);

    std::vector<occurrence_t> nucliobase_occurrence(reference.size());

    std::vector<std::thread> threads;

    for (int i = 0; i < num_of_threads; i++) {
        threads.push_back(std::thread(sequence_to_reference_map, std::ref(sequences),
                                      std::ref(reference), std::ref(minimizer_map), std::ref(minimizers),
                                      std::ref(minimizers_locations), std::ref(nucliobase_occurrence),
                                      kmer_size, minimizer_size, i));
    }

    for (auto &thread : threads) {
        thread.join();
    }

    std::ofstream rez_file;
    std::string rez_file_name("rez.report");
    rez_file.open(rez_file_name, std::fstream::out);

    for (int i = 0; i < reference.size(); i++) {
        occurrence_t occ = nucliobase_occurrence[i];
        //printf("%c  %d %d %d %d %d ", reference[i], occ.A, occ.C, occ.G, occ.T, occ.del);
        rez_file << reference[i]<<" " << occ.A << " " << occ.C << " " << occ.G << " " <<  occ.T << " " << occ.del;
        int j = 0;
        for (auto & a : occ.insert)
        {
        //    printf("  %d  %d %d %d %d  ", j, a.A, a.C, a.G, a.T);
            rez_file << "  " << j << "| " <<  a.A << " "  << a.C << " " <<  a.G << " " << a.T;
            j++;
        }
        //printf("\n");
        rez_file << std::endl;
    }
}
