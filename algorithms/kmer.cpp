#include "kmer.h"

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
                       std::vector<int> &minimizers,  std::vector<int> &minimizers_locations)
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
                minimizers_locations.push_back(i);
                unique_id++;
            }
            else
            {
                int id = it->second;
                minimizer_info_t info = {i, 1};
                index_map[id].push_back(info);
                minimizers.push_back(id);
                minimizers_locations.push_back(i);
            }
        }
        else
        {
            std::map<std::string, int>::iterator it =  minimizer_map.find(minimizer);
            int id = it->second;
            std::vector<minimizer_info_t> vec = index_map[id];
            vec.back().num_of_frames++;
        }

        prev_minimizer.clear();
        prev_minimizer.assign(minimizer);
    }
}

void minimize_sequence(std::string seq, int kmer_size, int minimizer_size,
                       std::map<std::string, int> &minimizer_map,
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

        if (minimizer.compare(prev_minimizer) != 0)
        {
            std::map<std::string, int>::iterator it =  minimizer_map.find(minimizer);
            if ( it ==  minimizer_map.end())
            {
                minimizers.push_back(unknown_minimizer_index);
                unknown_minimizer_index--;
            }
            else
            {
                int id = it->second;
                minimizers.push_back(id);
            }
        }
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

int read_reference(std::string& file_path, std::string& reference)
{
    std::ifstream file;
    char file_line[256];

    file.open(file_path, std::ifstream::in);

    if (!file.is_open()){
        printf("Error opening file %s",file_path.c_str());
        return -1;
    }

    file.getline(file_line, 256);

    while (true)
    {
        file.getline(file_line, 256);
        if (strlen(file_line) == 0)
            break;
        reference.append(file_line);
    }

    return 0;
}

int read_sequence(std::ifstream& sequences_file,std::string &sequence)
{
    static int sequence_line = 0;
    sequence_line++;
    char sequence_info[256] = {0};
    if (!sequences_file.is_open()){
        printf("File not open.");
        return -1;
    }

    sequences_file.getline(sequence_info, 256);

    if ( strlen(sequence_info) == 0)
    {
        sequences_file.close();
        return -2;
    }
    while(true)
    {
        char c = (char) sequences_file.get();
        if (c == '\n')
            break;
        sequence.push_back(c);
    }

    return sequence_line;
}


#ifdef KMER
int main()
#else
int kmer_main()
#endif
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
    std::string reference_file("data/ecoli.fasta");
    std::string sequence_file("data/ecoli_simulated_reads.fasta");

    int x = read_reference(reference_file, reference);

    if (x == -1)
        return 0;

    minimize_reference(reference, kmer_size, minimizer_size, minimizer_map, index_map, minimizers, minimizers_locations);

    sequences.open(sequence_file, std::ifstream::in);

    read_sequence(sequences, sequence);

    minimize_sequence(sequence, kmer_size, minimizer_size, minimizer_map, sequence_minimizers);

    int *a = &minimizers[0];
    int *b = &sequence_minimizers[0];

    subsequence_info_t rez = subsequence_size<int>(a, minimizers.size(), b, sequence_minimizers.size());

    printf("sequence original: len %ld minimizer indexes %ld %ld, original indexes %d %d \n",  rez.length, rez.first_index, rez.last_index, minimizers_locations[rez.first_index], minimizers_locations[rez.last_index]);

    create_complement(reverse, sequence);
    minimize_sequence(reverse, kmer_size, minimizer_size, minimizer_map, sequence_minimizers_reverse);

    a = &minimizers[0];
    b = &sequence_minimizers_reverse[0];

    rez = subsequence_size<int>(a, minimizers.size(), b, sequence_minimizers_reverse.size());

    printf("sequence reverse: len %ld minimizer indexes %ld %ld, original indexes %d %d \n",  rez.length, rez.first_index, rez.last_index, minimizers_locations[rez.first_index], minimizers_locations[rez.last_index]);
}
