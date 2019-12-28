#ifndef _HIRSCHBERG_H
#define _HIRSCHBERG_H

#include <string>
#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>

#define UP_LEFT 0
#define LEFT 1
#define UP 2

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


void needlemanWunschScore(const std::string& x, const std::string& y, std::vector<int>& upperLine, std::vector<int>& lowerLine, int sub, int del, int ins, int match);
std::pair<std::string, std::string>  needlemanWunsch(const std::string& x, const std::string& y, int ins, int del, int sub, int match);
std::pair<std::string, std::string> hirschbergAlgorithm(const std::string& x, const std::string& y, const std::string& reverseX, const std::string& reverseY, std::vector<int>& upperLine, std::vector<int>& lowerLine, std::vector<int>& bufferLine, int sub, int del, int ins, int match);
std::pair<std::string, std::string> hirschberg(const std::string& x, const std::string& y, int sub, int del, int ins, int match);

#endif //_HIRSCHBERG_H