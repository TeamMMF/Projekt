//
// Created by filip on 06.12.17..
//
#include "ProjectConfig.h"

#include <cstdio>
#include <memory>
#include "Common.hpp"
#include "FASTASampleClass.cpp"
#include "lcskpp.h"
#include "bioparser/bioparser.hpp"


using namespace std;

int main(int argc, char const *argv[]){
    const int k = 3;
    int length;

    //string sequence_file_path(argv[1]);
    string project_root(PROJECT_ROOT);

    vector<unique_ptr<FASTASampleClass>> fasta_reads;
    auto fasta_reader = bioparser::createReader<FASTASampleClass, bioparser::FastaReader>( project_root + "src/resources/sample.fasta");
    fasta_reader->read_objects(fasta_reads, static_cast<uint64_t>(-1));

    for (int i = 0, limit = fasta_reads.size(); i < limit; ++i) {
        for (int j = i+1; j < limit; ++j) {

            vector<pair<int,int>> rec;
            lcskpp_sparse_fast(fasta_reads[i]->get_data(),fasta_reads[j] -> get_data(), k, &length, &rec);

            printf("Length: %d\n",length);

            printf("Indices:");
            for(int k = 0, limit = rec.size(); k<limit; k++){
                printf(" (%d, %d)",rec[k].first, rec[k].second);
            }

            printf("\n");
        }
    }

    return 0;
}




