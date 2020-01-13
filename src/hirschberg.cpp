#include "hirschberg.h"


/*
    Called internally by the hirschberg_algorithm function.
    
    Calculates the Needleman-Wunsch alignment score for the given strings as part of the hirschberg algorithm.
*/
void needleman_wunsch_score(const std::string& x, const std::string& y, std::vector<int>& upper_line, std::vector<int>& lower_line, int sub, int del, int ins, int match) {
    upper_line[0] = 0;

    for(unsigned int i = 1, limit = y.size(); i <= limit; i++) {
        upper_line[i] = upper_line[i-1] + ins;
    }

    for(unsigned int i = 1, limit = x.size(); i <= limit; i++) {
        lower_line[0] = upper_line[0] + del;
        for(unsigned int j = 1, limit_2 = y.size(); j <= limit_2; j++) {
            if(x[i-1] == y[j-1]) {
                lower_line[j] = std::max({upper_line[j-1] + match, upper_line[j] + del, lower_line[j-1] + ins});
            } else {
                lower_line[j] = std::max({upper_line[j-1] + sub, upper_line[j] + del, lower_line[j-1] + ins});
            }
        }
        if(i <= limit - 1) {
            std::copy(lower_line.begin(), lower_line.begin() + y.size(), upper_line.begin());
        }
    }
}



/*
    Uses the Needleman-Wunsch algorithm to calculate global or half-global similarity.
    Used in the edge case when the length of one of the strings is 1.
    
    Called internally by the hirschberg_algorithm function.
    Can also be used externally.
*/
std::pair<std::string, std::string>  needleman_wunsch(const std::string& x, const std::string& y, int ins, int del, int sub, int match) {
    int height = x.size();
    int width = y.size();
    
    std::vector<std::vector<Cell>> cells(height + 1, std::vector<Cell>(width + 1));

    for(int i = 0; i <= width; i++) {
        cells[0][i].cost = -i;
        cells[0][i].parent = LEFT;
    }

    for(int i = 0; i <= height; i++) {
        cells[i][0].cost = -i;
        cells[i][0].parent = UP;
    }

    for(unsigned int i = 1; i <= x.size(); i++) {
        for(unsigned int j = 1; j <= y.size(); j++) {

            int match_cost;
            if(x[i-1] == y[j-1]) {
                match_cost = cells[i-1][j-1].cost + match;
            } else {
                match_cost = cells[i-1][j-1].cost + sub;
            }
            int insCost = cells[i][j-1].cost + ins;
            int delCost = cells[i-1][j].cost + del;
            
            update_cell(cells[i][j], match_cost, insCost, delCost, true);
        }
    }

    std::string z;
    std::string w;
    unsigned int i = x.size();
    unsigned int j = y.size();
    Cell currentCell = cells[i][j];
    while(i != 0 || j != 0) {
        if(currentCell.parent == UP_LEFT) {
            z.append(1,x[i-1]);
            w.append(1,y[j-1]);
            i--;
            j--;
        } else if (currentCell.parent == UP) {
            z.append(1,x[i-1]);
            w.append(1,'-');
            i--;
        } else {
            z.append(1,'-');
            w.append(1,y[j-1]);
            j--;
        }
        currentCell = cells[i][j];
    }

    std::reverse(z.begin(), z.end());
    std::reverse(w.begin(), w.end());
    return std::make_pair(z, w);
}

/*
    Internally called by the hirschberg function.
    
    Implements the recursive portion of the algorithm.
*/
std::pair<std::string, std::string> hirschberg_algorithm(const std::string& x, const std::string& y, const std::string& reverseX, const std::string& reverseY, std::vector<int>& upper_line, std::vector<int>& lower_line, std::vector<int>& buffer_line, int sub, int del, int ins, int match) {
    auto xlen = x.size();
    auto xmid = xlen / 2;
    auto ylen = y.size();

    std::string z;
    std::string w;

    if(xlen == 0) {
        for(unsigned int i = 0; i < ylen; i++) {
            z.append(1,'-');
            w.append(1,y[i]);
        }
        return std::make_pair(z, w);
    }
    if(ylen == 0) {
        for(unsigned int i = 0; i < xlen; i++) {
            z.append(1,x[i]);
            w.append(1,'-');
        }
        return std::make_pair(z, w);
    }
    if(xlen == 1 || ylen == 1) {
        return needleman_wunsch(x, y, sub, del, ins, match);
    }

    needleman_wunsch_score(x.substr(0, xmid), y, upper_line, lower_line, sub, del, ins, match);
    std::copy(lower_line.begin(), lower_line.begin() + ylen + 1, buffer_line.begin());
    int reversemid = xlen - xmid;
    needleman_wunsch_score(reverseX.substr(0, reversemid), reverseY, upper_line, lower_line, sub, del, ins, match);
    

    int argmax = 0;
    int max = buffer_line[0] + lower_line[ylen];
    for(unsigned int i = 1; i <= ylen; i++) {
        if(buffer_line[i] + lower_line[ylen - i] > max) {
            max = buffer_line[i] + lower_line[ylen - i];
            argmax = i;
        }
    }

    int reverseargmax = ylen - argmax;
    auto left = hirschberg_algorithm(x.substr(0, xmid), y.substr(0, argmax), reverseX.substr(reversemid, xlen), reverseY.substr(reverseargmax, ylen), upper_line, lower_line, buffer_line, sub, del, ins, match);
    auto right = hirschberg_algorithm(x.substr(xmid, xlen), y.substr(argmax, ylen), reverseX.substr(0, reversemid), reverseY.substr(0, reverseargmax), upper_line, lower_line, buffer_line, sub, del, ins, match);
    return std::make_pair(left.first + right.first, left.second + right.second);
}

/*
    Aligns the x and y strings using the Hirschberg algorithm with the given costs.
*/
std::pair<std::string, std::string> hirschberg(const std::string& x, const std::string& y, int sub, int del, int ins, int match) {
    std::vector<int> upper_line(y.size() + 1);
    std::vector<int> lower_line(y.size() + 1);
    std::vector<int> buffer_line(y.size() + 1);

    std::string reverseX;
    reverseX.assign(x.rbegin(), x.rend());
    std::string reverseY;
    reverseY.assign(y.rbegin(), y.rend());
    return hirschberg_algorithm(x, y, reverseX, reverseY, upper_line, lower_line, buffer_line, sub, del, ins, match);
}

#ifdef HIRSCHBERG
int main(void) {
#else 
int hirsch_main(int argc, char *argv[]) {
#endif
    std::pair<std::string, std::string> r = hirschberg("AGTAACGCA", "TATGC", -1, -2, -2, 2);
    std::cout << r.first << ',' << r.second << '\n';
    r = needleman_wunsch("AGTAACGCA", "TATGC", -1, -1, -2, 2);
    std::cout << r.first << ',' << r.second << '\n';
    return 0;
}