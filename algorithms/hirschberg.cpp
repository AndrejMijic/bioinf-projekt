#include <string>
#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>

#define UP_LEFT 0
#define LEFT 1
#define UP 2

void needlemanWunschScore(const std::string& x, const std::string& y, std::vector<int>& upperLine, std::vector<int>& lowerLine, int sub, int del, int ins, int match) {
	upperLine[0] = 0;
	//std::cout << upperLine[0] << ' ';
	for(unsigned int i = 1, limit = y.size(); i <= limit; i++) {
		upperLine[i] = upperLine[i-1] + ins;
		//std::cout << upperLine[i] << ' ';
	}
	//std::cout << '\n';
	//for(unsigned int i = 0; )
	for(unsigned int i = 1, limit = x.size(); i <= limit; i++) {
		lowerLine[0] = upperLine[0] + del;
		//std::cout << lowerLine[0] << ' ';
		for(unsigned int j = 1, limit_2 = y.size(); j <= limit_2; j++) {
			if(x[i-1] == y[j-1]) {
				lowerLine[j] = std::max({upperLine[j-1] + match, upperLine[j] + del, lowerLine[j-1] + ins});
			} else {
				lowerLine[j] = std::max({upperLine[j-1] + sub, upperLine[j] + del, lowerLine[j-1] + ins});
			}
			//std::cout << lowerLine[j] << ' ';
		}
		//std::cout << '\n';
		if(i <= limit - 1) {
			std::copy(lowerLine.begin(), lowerLine.begin() + y.size(), upperLine.begin());
		}
	}
}


/*
	Dynamic programming table cell.
*/
struct Cell {
	int cost;
	int parent;

	Cell(int cost, int parent) : cost{cost}, parent{parent} {}
	Cell() : cost{0}, parent{0} {}
};

/*
	Sets the value of the cell depending on the costs for matching, insertion and deletion.
*/
inline void updateCell(Cell& cell, int match, int ins, int del, bool similarity = false) {

	cell.cost = match;
	cell.parent = UP_LEFT;

	if(similarity) {
		if(ins > cell.cost) {
			cell.cost = ins;
			cell.parent = LEFT;
		} 
		if(del > cell.cost) {
			cell.cost = del;
			cell.parent = UP;
		}
	} else {
		if(ins < cell.cost) {
			cell.cost = ins;
			cell.parent = LEFT;
		} 
		if(del < cell.cost) {
			cell.cost = del;
			cell.parent = UP;
		}	
	}
}


/*
	Uses the Needleman-Wunsch algorithm to calculate global or half-global similarity.
*/
std::pair<std::string, std::string>  needlemanWunsch(const std::string& x, const std::string& y, int ins, int del, int sub, int match) {
	int height = x.size();
	int width = y.size();

	//std::cout << height;
	//std::cout << width;
	
	std::vector<std::vector<Cell>> cells(height + 1, std::vector<Cell>(width + 1));

	//std::cout << cells.size();
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

			//int match = cells[i-1][j-1].cost + costs[decode(s[i-1])][decode(t[j-1])];
			int matchCost;
			if(x[i-1] == y[j-1]) {
				matchCost = cells[i-1][j-1].cost + match;
			} else {
				matchCost = cells[i-1][j-1].cost + sub;
			}
			int insCost = cells[i][j-1].cost + ins;
			int delCost = cells[i-1][j].cost + del;
			
			updateCell(cells[i][j], matchCost, insCost, delCost, true);

			//std::cout << "Putting: " << cells[i][j].cost << ", on (" << i << ", " << j <<"), with direction: " << cells[i][j].parent << "\n";
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


std::pair<std::string, std::string> hirschbergAlgorithm(const std::string& x, const std::string& y, const std::string& reverseX, const std::string& reverseY, std::vector<int>& upperLine, std::vector<int>& lowerLine, std::vector<int>& bufferLine, int sub, int del, int ins, int match) {
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
		return needlemanWunsch(x, y, sub, del, ins, match);
	}

	needlemanWunschScore(x.substr(0, xmid), y, upperLine, lowerLine, sub, del, ins, match);
	std::copy(lowerLine.begin(), lowerLine.begin() + ylen + 1, bufferLine.begin());
	int reversemid = xlen - xmid;
	/*if(xlen % 2 == 0) {
		needlemanWunschScore(reverseX.substr(0, xmid), reverseY, upperLine, lowerLine, sub, del, ins, match);
	} else {
		needlemanWunschScore(reverseX.substr(0, xmid + 1), reverseY, upperLine, lowerLine, sub, del, ins, match);
	}*/
	needlemanWunschScore(reverseX.substr(0, reversemid), reverseY, upperLine, lowerLine, sub, del, ins, match);
	

	int argmax = 0;
	int max = bufferLine[0] + lowerLine[ylen];
	/*for(unsigned int i = 0; i <= ylen; i++) {
		std::cout << bufferLine[i] << ' ';
	}
	std::cout << '\n';
	for(unsigned int i = 0; i <= ylen; i++) {
		std::cout << lowerLine[ylen - i] << ' ';
	}
	std::cout << '\n';*/
	for(unsigned int i = 1; i <= ylen; i++) {
		if(bufferLine[i] + lowerLine[ylen - i] > max) {
			max = bufferLine[i] + lowerLine[ylen - i];
			argmax = i;
		}
	}

	int reverseargmax = ylen - argmax;
	auto left = hirschbergAlgorithm(x.substr(0, xmid), y.substr(0, argmax), reverseX.substr(reversemid, xlen), reverseY.substr(reverseargmax, ylen), upperLine, lowerLine, bufferLine, sub, del, ins, match);
	auto right = hirschbergAlgorithm(x.substr(xmid, xlen), y.substr(argmax, ylen), reverseX.substr(0, reversemid), reverseY.substr(0, reverseargmax), upperLine, lowerLine, bufferLine, sub, del, ins, match);
	return std::make_pair(left.first + right.first, left.second + right.second);
}


std::pair<std::string, std::string> hirschberg(const std::string& x, const std::string& y, int sub, int del, int ins, int match) {
	std::vector<int> upperLine(y.size() + 1);
	std::vector<int> lowerLine(y.size() + 1);
	std::vector<int> bufferLine(y.size() + 1);

	std::string reverseX;
	reverseX.assign(x.rbegin(), x.rend());
	std::string reverseY;
	reverseY.assign(y.rbegin(), y.rend());
	return hirschbergAlgorithm(x, y, reverseX, reverseY, upperLine, lowerLine, bufferLine, sub, del, ins, match);
}


int main(int argc, char *argv[]) {
	std::pair<std::string, std::string> r = hirschberg("AGTAACGCA", "TATGC", -1, -2, -2, 2);
	std::cout << r.first << ',' << r.second << '\n';
	r = needlemanWunsch("AGTAACGCA", "TATGC", -1, -1, -2, 2);
	std::cout << r.first << ',' << r.second << '\n';
	return 0;
}