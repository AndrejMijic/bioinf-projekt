#include <cstdlib>
#include <string>
#include <ctime>
#include <fstream>
#include <iostream>

int main(int argc, char* argv[]) {

	/*if(argc != 3) {
		return 1;
	} */

	int length = atoi(argv[1]);
	std::string reference;
	std::srand(std::time(nullptr));

	for(int i = 0; i < length; i++) {
		int value = std::rand() % 4;
		switch(value) {
			case 0:
				reference.push_back('A');
				break;
			case 1:
				reference.push_back('C');
				break;
			case 2:
				reference.push_back('T');
				break;
			case 3:
				reference.push_back('G');
				break;
			default:
				break;											
		}
	}

	for(int i = 0; i < reference.size(); i++) {
		int value = std::rand() % 20;
		if(value == 0) {
			value = std::rand() % 9;
			try{
			switch(value) {
				case 0:
					reference[i] ='A';
					break;
				case 1:
					reference[i] ='C';
					break;
				case 2:
					reference[i] ='T';
					break;
				case 3:
					reference[i] ='G';
					break;
				case 4:
				value = std::rand() % 100;
					if(value == 0) {
					reference.insert(i, "A");
				}
					break;
				case 5:
				value = std::rand() % 100;
					if(value == 0) {
					reference.insert(i, "C");
				}
					break;
				case 6:
				value = std::rand() % 100;
					if(value == 0) {
					reference.insert(i, "G");
				}
					break;
				case 7:
				value = std::rand() % 100;
					if(value == 0) {
					reference.insert(i, "T");
				}
					break;
				case 8:
					value = std::rand() % 100;
					if(value == 0) {
						reference.erase(i);
					}
				default:
					break;											
			}
			} catch	(const std::out_of_range& e) {
				std::cout << "Error at index " << i << '\n';
				std::cout << "Value was " << value << '\n';
				std::cout << "String length was " << reference.size() << '\n';
				return 1;
			}
		}
	}

	length = reference.size();

	std::ofstream os(argv[2], std::ios::out);
	os << ">test\n";
	os << reference << '\n';

	std::ofstream os_read(argv[3], std::ios::out);
	for(int i = 0; i < 800;) {
		int start = std::rand() % length;
		int end = std::rand() % length;

		if(start > end) {
			int tmp = start;
			start = end;
			end = tmp;
		}

		//std::cout << start << '\n';
		//std::cout << end << "\n\n";
		if(end - start < length / 50 || end - start > length / 20) {
			continue;
		}

		std::string read;
		for(int j = start; j < end; j++) {
			read.push_back(reference[j]);
		}

		/*for(int j = 0; j < read.size(); j++) {
			int value = std::rand() % 20;
			if(value == 0) {
				value = std::rand() % 4;
				switch(value) {
					case 0:
						read[j] ='A';
						break;
					case 1:
						read[j] ='C';
						break;
					case 2:
						read[j] ='T';
						break;
					case 3:
						read[j] ='G';
						break;
					default:
						break;											
				}	
			}		
		}*/

		os_read << ">simulated " << start << "to" << end <<'\n';
		os_read << read << '\n';
		i++;

		
	} 

	return 0;
}