//
// Created by matep on 07/11/2017.
//

#include "ProjectConfig.h"

#include <string>
#include <vector>
#include <iostream>
#include <unordered_set>
#include <memory>
#include <chrono>
#include <iomanip>

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
        //vector<triplet> vector = find_minimizers(4, 3, "ACGACTGGTCAGAGT");

        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
        vector<triplet> vector1 = find_minimizers(5, 15, (*fasta_object).get_data());
        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        /*for (auto &i : vector) {
            triplet t = static_cast<tuple<string, int, string> &&>(i);
            cout << get<0>(t) << " " << get<1>(t) << " " << get<2>(t) << endl;
        }*/

        //auto vector2 = find_minimizers2(4, 3, "ACGACTGGTCAGAGT");
        chrono::high_resolution_clock::time_point t3 = chrono::high_resolution_clock::now();
        auto vector2 = find_minimizers2(5, 15, (*fasta_object).get_data());
        chrono::high_resolution_clock::time_point t4 = chrono::high_resolution_clock::now();
        cout << vector2.size() << endl;
        for (auto &i : vector2) {
            //triple t = static_cast<tuple<uint64_t , int, int> &&>(i);
            std::cout <<setw(10) << std::get<0>(i) << " " << setw(offset_width)<< std::get<1>(i) << " " << std::get<2>(i) << endl;
        }

        auto duration1 = chrono::duration_cast<chrono::microseconds>( t2 - t1 ).count();
        auto duration2 = chrono::duration_cast<chrono::microseconds>( t4 - t3 ).count();

        cout << duration1 << endl;
        cout << duration2 << endl;
        break;
    }

    return 0;
}

