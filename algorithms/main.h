typedef struct occurrence {
    int A = 0;
    int C = 0;
    int G = 0;
    int T = 0;
    int del = 0;
    std::vector<struct occurrence> insert;
} occurrence_t;