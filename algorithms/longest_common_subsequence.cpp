#include "longest_common_subsequence.h"

//for testing of correctness  http://lcs-demo.sourceforge.net/
int main(void)
{
    char a[100] = "FEHAFGJIGCAJIGBHGEBJHECHAFIJHHIGFIBFEAEFIDGJIAECHDIIHEAHGFICDIDFGCACCCDBADAGAHAECJBCFDAAAHHDJGGGIAEC";
    char b[100] = "DIIHEAHHGFICIDFGCACCCGDBADAGCAHAEJBCFDAAAHHDJG";

    char c[77] = "DADCDCDADBB";
    char d[33] = "DBCADCACABD";
    subsequence_info_t rez = subsequence_size<char>(a, 100, b, 46);
    printf("%ld %ld %ld\n", rez.length, rez.first_index, rez.last_index);
    rez = subsequence_size<char>(c, 11, d, 11);
    printf("%ld %ld %ld\n",  rez.length, rez.first_index, rez.last_index);
}

