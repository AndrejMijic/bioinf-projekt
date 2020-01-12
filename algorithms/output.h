#ifndef _OUTPUT_H
#define _OUTPUT_H
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#include "main.h"

void write_to_CSV(std::string filename, std::vector<occurrence_t>& occurences, int minCoverage, std::string& refGenome);
#endif

