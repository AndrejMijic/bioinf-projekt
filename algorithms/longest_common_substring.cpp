#include "longest_common_substring.h"



int main(void)
{
    char a[1002] = "AACCDCBBDDBBABBCDDCB";
    char b[1002] = "BBABBACDCCB";

    substring_info_t rez = substring_size<char>(&a[0], 20, &b[0], 8);
    printf("%ld %ld \n", rez.length, rez.first_index);

}

