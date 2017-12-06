//
// Created by filip on 06.12.17..
//

#include "ProjectConfig.h"

#include <iostream>
#include <memory>
#include "Common.hpp"
#include "FASTASampleClass.cpp"
#include "lcskpp.h"
#include "bioparser/bioparser.hpp"

#define K 3
#define W 5

vector<pair<int,int>> get_matches(vector<tuple<uint64_t, int, int>> minimizers1, vector<tuple<uint64_t, int, int>> minimizers2){

    vector<pair<int,int>> matches;
    for(int i= 0, li = minimizers1.size(); i<li; i++){
        tuple<uint64_t, int, int> min1 = minimizers1[i];
        for(int j = 0, lj = minimizers2.size(); j<lj; j++){
            tuple<uint64_t, int, int> min2 = minimizers2[j];
            if(get<2>(min2)==1 || get<2>(min1)==1){
                continue;
            }
            if(get<0>(min2) == get<0>(min1)){
                matches.push_back(make_pair(i,j));
            }
        }
    }
    return matches;
};

int main(int argc, char const *argv[]){
    int length;

    //string sequence_file_path(argv[1]);
    string project_root(PROJECT_ROOT);

    vector<unique_ptr<FASTASampleClass>> fasta_reads;
    auto fasta_reader = bioparser::createReader<FASTASampleClass, bioparser::FastaReader>( project_root + "src/resources/sample.fasta");
    fasta_reader->read_objects(fasta_reads, static_cast<uint64_t>(-1));

    vector<vector<tuple<uint64_t, int, int>>> minimizers_for_read;

    for (int i = 0, limit = fasta_reads.size(); i < limit; ++i) {
        vector<tuple<uint64_t, int, int>> minimizers = find_minimizers2(W, K, fasta_reads[i] -> get_data());
        minimizers_for_read.push_back(minimizers);
    }

    for (int i = 0, limit = fasta_reads.size(); i < limit; ++i) {
        for (int j = i+1; j < limit; ++j) {
            vector<pair<int,int>> mathces = get_matches(minimizers_for_read[i],minimizers_for_read[j]);

            vector<pair<int,int>> rec;
            lcskpp_sparse_fast(mathces, K, &length, &rec);
            cout <<length<<endl;
            string data1 = fasta_reads[i] -> get_data();
            string data2 = fasta_reads[j] -> get_data();
            for(int k = 0, limit = rec.size(); k<limit; k++){
                cout <<"(" << rec[k].first<<","<< rec[k].second<<")";
            }
            cout <<endl;
        }
    }

    return 0;
}

