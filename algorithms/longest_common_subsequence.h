#ifndef _SUBSEQUENCE_H
#define _SUBSEQUENCE_H
#include <malloc.h>
#include <stdio.h>

typedef struct subsequence_info
{
    int length;
    int first_index;
    int last_index;
} subsequence_info_t;

subsequence_info_t *subsequence_info_max(subsequence_info_t *a, subsequence_info_t *b)
{
    if (a->length >= b->length)
        return a;
    else
        return b;
}


subsequence_info_t *subsequence_info_max(subsequence_info_t *a, subsequence_info_t *b);

//because of template
template<typename T>
subsequence_info_t subsequence_size(T *ref, int ref_size, T *seq, int seq_size)
{
    subsequence_info_t *ar[2];
    long i, j;
    ar[0] = (subsequence_info_t *)malloc((ref_size + 1) * sizeof(subsequence_info_t));
    ar[1] = (subsequence_info_t *)malloc((ref_size + 1)* sizeof(subsequence_info_t));
    
    for (i = 0; i< ref_size + 1; i++)
    {
        ar[0][i].length = 0;
        ar[0][i].first_index = -1;
        ar[0][i].last_index = -1;
    }
    
    ar[1][0].length = 0;
    ar[1][0].first_index = -1;
    ar[1][0].last_index = -1;

    subsequence_info_t ret;
    int top_index = 1;
    int bottom_index = 0;
    subsequence_info_t *helper;

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
                    ar[bottom_index][j + 1].last_index = j + 1;
                    ar[bottom_index][j + 1].length = 1;
                }
                else
                {
                    ar[bottom_index][j + 1].length = 1 + ar[top_index][j].length;
                    ar[bottom_index][j + 1].last_index = j + 1;
                    ar[bottom_index][j + 1].first_index = ar[top_index][j].first_index;
                    int x = 9;
                }
                
            }
            else
            {
                helper = subsequence_info_max(&(ar[top_index][j + 1]), &(ar[bottom_index][j]));
                ar[bottom_index][j + 1].length = helper->length;
                ar[bottom_index][j + 1].first_index = helper->first_index;
                ar[bottom_index][j + 1].last_index = helper->last_index;
                int x = 9;
            }
        }
    }

    ret = ar[bottom_index][ref_size];
    ret.first_index--; // index 0 is for empty array
    ret.last_index--;
    free(ar[0]);
    free(ar[1]);

    return ret;
}

#endif // _SUBSEQUENCE_H