#include "longest_common_subsequence.h"

//for testing of correctness  http://lcs-demo.sourceforge.net/
int main(void)
{
    char a[121] = "FEHAFGJIGCAJIGBHGEBJHECHAFIJHHIGFIBFEAEFIDGJIAECHDIIHEAHGFICDIDFGCACCCDBADAGAHAECJBCFDAAAHHDJGGGIAEC";
    char b[100] = "DIIHEAHHGFICIDFGCACCCGDBADAGCAHAEJBCFDAAAHHDJG";

    char c[77] = "DADCDCDADBB";
    char d[33] = "DBCADCACABD";

    subsequence_info_t rez = subsequence_size<char>( b, 46,a, 120);
    printf("%d %d %d %d %d\n", rez.length, rez.first_index_ref, rez.last_index_ref, rez.first_index_seq, rez.last_index_seq);
    rez = subsequence_size<char>(c, 11, d, 11);
    printf("%d %d %d %d %d\n", rez.length, rez.first_index_ref, rez.last_index_ref, rez.first_index_seq, rez.last_index_seq);
}

