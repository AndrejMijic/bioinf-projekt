#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

#include "main.h"

void write_to_CSV(std::string filename, std::vector<occurrence_t>& occurences, int min_coverage, std::string& ref_genome);

#endif // OUTPUT_H_

