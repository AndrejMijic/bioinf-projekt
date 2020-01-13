#ifndef HIRSCHBERG_H_
#define HIRSCHBERG_H_

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
inline void update_cell(Cell& cell, int match, int ins, int del, bool similarity = false) {

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


void needleman_wunsch_score(const std::string& x, const std::string& y, std::vector<int>& upperLine, std::vector<int>& lowerLine, int sub, int del, int ins, int match);
std::pair<std::string, std::string>  needleman_wunsch(const std::string& x, const std::string& y, int ins, int del, int sub, int match);
std::pair<std::string, std::string> hirschberg_algorithm(const std::string& x, const std::string& y, const std::string& reverseX, const std::string& reverseY, std::vector<int>& upperLine, std::vector<int>& lowerLine, std::vector<int>& bufferLine, int sub, int del, int ins, int match);
std::pair<std::string, std::string> hirschberg(const std::string& x, const std::string& y, int sub, int del, int ins, int match);

#endif // HIRSCHBERG_H_