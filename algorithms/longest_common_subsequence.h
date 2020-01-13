#ifndef SUBSEQUENCE_H_
#define SUBSEQUENCE_H_
#include <malloc.h>
#include <stdio.h>


/*
    Structure that holds the longest common subsequence result.
*/
struct subsequence_info_t
{
    int length;
    int first_index_ref;
    int last_index_ref;
    int first_index_seq;
    int last_index_seq;
};

/*
    Calculates the longest common subsequence between the given strings.
*/
template<typename T>
subsequence_info_t subsequence_size(T *ref, int ref_size, T *seq, int seq_size)
{
    subsequence_info_t *ar[2];
    long i, j;
    int seq_index = 0;
    ar[0] = (subsequence_info_t *)malloc((ref_size + 1) * sizeof(subsequence_info_t));
    ar[1] = (subsequence_info_t *)malloc((ref_size + 1)* sizeof(subsequence_info_t));

    for (i = 0; i< ref_size + 1; i++)
    {
        ar[0][i].length = 0;
        ar[0][i].first_index_ref = -1;
        ar[0][i].last_index_ref = -1;
        ar[0][i].first_index_seq = -1;
        ar[0][i].last_index_seq = -1;
    }

    ar[1][0].length = 0;
    ar[1][0].first_index_ref = -1;
    ar[1][0].last_index_ref = -1;
    ar[1][0].first_index_seq = -1;
    ar[1][0].last_index_seq = -1;

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
                    ar[bottom_index][j + 1].first_index_ref = j + 1;
                    ar[bottom_index][j + 1].last_index_ref = j + 1;
                    ar[bottom_index][j + 1].first_index_seq = seq_index;
                    ar[bottom_index][j + 1].last_index_seq = seq_index;
                    ar[bottom_index][j + 1].length = 1;
                }
                else
                {
                    ar[bottom_index][j + 1].length = 1 + ar[top_index][j].length;
                    ar[bottom_index][j + 1].last_index_ref = j + 1;
                    ar[bottom_index][j + 1].first_index_ref = ar[top_index][j].first_index_ref;
                    ar[bottom_index][j + 1].first_index_seq = ar[top_index][j].first_index_seq;
                    ar[bottom_index][j + 1].last_index_seq = seq_index;
                }
            }
            else
            {
                helper = ar[top_index][j + 1].length >= ar[bottom_index][j].length ? &(ar[top_index][j + 1]) : &(ar[bottom_index][j]);
                ar[bottom_index][j + 1].length = helper->length;
                ar[bottom_index][j + 1].first_index_ref = helper->first_index_ref;
                ar[bottom_index][j + 1].last_index_ref = helper->last_index_ref;
                ar[bottom_index][j + 1].first_index_seq = helper->first_index_seq;
                ar[bottom_index][j + 1].last_index_seq = helper->last_index_seq;
            }
        }
        seq_index++;
    }

    ret = ar[bottom_index][ref_size];
    ret.first_index_ref--; // index 0 is for empty array
    ret.last_index_ref--;
    free(ar[0]);
    free(ar[1]);

    return ret;
}

#endif // SUBSEQUENCE_H_