//
// Created by matep on 07/11/2017.
//

#include "ProjectConfig.h"

#include <string>
#include <vector>
#include <iostream>
#include <unordered_set>
#include <memory>

#include "FASTASampleClass.cpp"
#include "bioparser/bioparser.hpp"
#include "Common.hpp"

using namespace std;


int main(int argc, char const *argv[]) {

    string project_root(PROJECT_ROOT);

    vector<unique_ptr<FASTASampleClass>> fasta_objects;
    auto fasta_reader = bioparser::createReader<FASTASampleClass, bioparser::FastaReader>(
            project_root + "src/resources/sample.fasta");

    fasta_reader->read_objects(fasta_objects, static_cast<uint64_t>(-1));

    for (auto &fasta_object : fasta_objects) {

        //vector<string> vector = find_kmer(100, (*fasta_object).get_data());
        vector<triplet> vector = find_minimizers(4, 3, "123456354748");
        for (auto &i : vector) {
            triplet t = static_cast<tuple<string, int, string> &&>(i);
            cout << get<0>(t) << " " << get<1>(t) << " " << get<2>(t) << endl;
        }
        break;
    }

    return 0;
}

