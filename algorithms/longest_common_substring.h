#ifndef _SUBSTRING_H
#define _SUBSTRING_H
#include <malloc.h>
#include <stdio.h>


typedef struct substring_info
{
    long length;
    long first_index;
} substring_info_t;



//because of template
template<typename T>
substring_info_t substring_size( T *ref, long ref_size, T *seq, long seq_size)
{
    substring_info_t *ar[2];
    long i, j;
    ar[0] = (substring_info_t *)malloc((ref_size + 1) * sizeof(substring_info_t));
    ar[1] = (substring_info_t *)malloc((ref_size + 1)* sizeof(substring_info_t));
    
    for (i = 0; i< ref_size + 1; i++)
    {
        ar[0][i].length = 0;
        ar[0][i].first_index = -1;
    }
    
    ar[1][0].length = 0;
    ar[1][0].first_index = -1;

    substring_info_t ret;
    int top_index = 1;
    int bottom_index = 0;
    substring_info_t helper = {0,0};

    for (i = 0; i < seq_size; i++)
    {
        top_index = bottom_index;
        bottom_index = top_index == 0 ? 1 : 0;

        for (j = 0; j < ref_size; j++)
        {
            if (ref[j] == seq[i])
            {
                if (ar[top_index][j].length == 0)
                {
                    ar[bottom_index][j + 1].first_index = j + 1;
                    ar[bottom_index][j + 1].length = 1;
                    if ( 1 > helper.length){
                        helper.first_index = j + 1;
                        helper.length = 1;
                    }
                }
                else
                {
                    ar[bottom_index][j + 1].length = 1 + ar[top_index][j].length;
                    ar[bottom_index][j + 1].first_index = ar[top_index][j].first_index;
                }

                if ( ar[bottom_index][j + 1].length > helper.length)
                {
                    helper.first_index = ar[bottom_index][j + 1].first_index;
                    helper.length = ar[bottom_index][j + 1].length;
                }
                
            }
            else
            {
                ar[bottom_index][j + 1].length = 0;
                ar[bottom_index][j + 1].first_index = -1;
            }
        }
    }

    ret = helper;
    ret.first_index--; // index 0 is for empty array
    free(ar[0]);
    free(ar[1]);

    return ret;
}

#endif //_SUBSTRING_H