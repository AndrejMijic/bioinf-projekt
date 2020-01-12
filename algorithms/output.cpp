#include "output.h"


//using namespace std;



void write_to_CSV(std::string filename, std::vector<occurrence_t>& occurences, int minCoverage, std::string& refGenome){
    std::ofstream outputFile;
    outputFile.open(filename, std::fstream::out);

    for (int i=0; i < occurences.size(); i++){
        
        if (occurences[i].A == 0 && occurences[i].C == 0 && occurences[i].G == 0 && occurences[i].T == 0 && occurences[i].del == 0) {continue;}

        if (occurences[i].A > occurences[i].C && occurences[i].A > occurences[i].G  && occurences[i].A > occurences[i].T && occurences[i].A > occurences[i].del){
            if (refGenome[i] != 'A' && occurences[i].A > minCoverage){
                outputFile << 'X' << ',' << i << ',' << 'A' << std::endl;
            }
        } else if (occurences[i].C > occurences[i].A && occurences[i].C > occurences[i].G  && occurences[i].C > occurences[i].T && occurences[i].C > occurences[i].del){
            if (refGenome[i] != 'C' && occurences[i].C > minCoverage){
                outputFile << 'X' << ',' << i << ',' << 'C' << std::endl;
            }
        } else if (occurences[i].G > occurences[i].A && occurences[i].G > occurences[i].C  && occurences[i].G > occurences[i].T && occurences[i].G > occurences[i].del){
            if (refGenome[i] != 'G' && occurences[i].G > minCoverage){
                outputFile << 'X' << ',' << i << ',' << 'G' << std::endl;
            }
        } else if (occurences[i].T > occurences[i].A && occurences[i].T > occurences[i].C  && occurences[i].T > occurences[i].G && occurences[i].T > occurences[i].del){
            if (refGenome[i] != 'T' && occurences[i].T > minCoverage){
                outputFile << 'X' << ',' << i << ',' << 'T' << std::endl;
            }
        } else if (occurences[i].del > occurences[i].A && occurences[i].del > occurences[i].C  && occurences[i].del > occurences[i].G && occurences[i].del > occurences[i].T){
            if(occurences[i].del > minCoverage)
                outputFile << 'D' << ',' << i << ',' << '-' << std::endl;
            
        }

        for (occurrence_t& ins: occurences[i].insert){
            if (ins.A > ins.C && ins.A > ins.G  && ins.A > ins.T){
                if (ins.A > minCoverage)
                    outputFile << 'I' << ',' << i + 1 << ',' << 'A' << std::endl;

            } else if (ins.C > ins.A && ins.C > ins.G  && ins.C > ins.T){
                if (ins.C > minCoverage)
                    outputFile << 'I' << ',' << i + 1 << ',' << 'C' << std::endl;

            } else if (ins.G > ins.A && ins.G > ins.C  && ins.G > ins.T){
                if (ins.G > minCoverage)
                    outputFile << 'I' << ',' << i + 1 << ',' << 'G' << std::endl;

            } else if (ins.T > ins.A && ins.T > ins.C  && ins.T > ins.G){
                if (ins.T > minCoverage)
                    outputFile << 'I' << ',' << i + 1 << ',' << 'T' << std::endl;
                
            }     
        }
        
    }

    outputFile.close();

}

