//
// Created by matep on 07/11/2017.
//

#include "ProjectConfig.h"

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>
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
            project_root + "src/resources/lambda_reads.fasta");

    fasta_reader->read_objects(fasta_objects, static_cast<uint64_t>(-1));

    /*
    vector<unique_ptr<FASTASampleClass>> fasta_objects2;
    auto fasta_reader2 = bioparser::createReader<FASTASampleClass, bioparser::FastaReader>(
            project_root + "src/resources/lambda_reference.fasta");

    fasta_reader2->read_objects(fasta_objects2, static_cast<uint64_t>(-1));

    for (auto &fasta_object : fasta_objects) {

        //vector<string> vector = find_kmer(100, (*fasta_object).get_data());
        //vector<triplet> vector = find_minimizers(4, 3, "ACGACTGGTCAGAGT");

        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
        vector<triplet> vector1 = find_minimizers(5, 15, (*fasta_object).get_data());
        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        /for (auto &i : vector) {
         //   triplet t = static_cast<tuple<string, int, string> &&>(i);
        //    cout << get<0>(t) << " " << get<1>(t) << " " << get<2>(t) << endl;
       // }

        //auto vector2 = find_minimizers2(4, 3, "ACGACTGGTCAGAGT");
        chrono::high_resolution_clock::time_point t3 = chrono::high_resolution_clock::now();
        auto vector2 = find_minimizers2(5, 15, (*fasta_object).get_data());
        chrono::high_resolution_clock::time_point t4 = chrono::high_resolution_clock::now();
        cout << vector2.size() << endl;
        uint64_t  maks1 = 0;
        for (auto &i : vector2) {
            maks1 = maks1 < get<0>(i) ? get<0>(i) : maks1;
            std::cout <<setw(10) << std::get<0>(i) << " " << setw(offset_width)<< std::get<1>(i) << " " << std::get<2>(i) << endl;
        }

        vector<string> v;
        for(auto &f : fasta_objects){
            string s = (*f).get_data();
            v.push_back(s);
        }

        /*
        //vector<string> v;
        //string s = (*fasta_object).get_data();
        //v.push_back(s);


        unordered_multimap<uint64_t, tuple<string, int, int>, function<size_t( uint64_t)>> multiset = indexSequences(v, 5, 15);

        uint64_t  maks2 = 0;
       /*
        //for (auto &i : multiset) {
        //   maks2 = maks2 < get<0>(i) ? get<0>(i) : maks2;
        //    std::cout <<setw(10) << i.first << " " << setw(offset_width)<< std::get<0>(i.second) << " " << setw(offset_width) << std::get<1>(i.second) << " " << std::get<2>(i.second) << endl;
        //}

        cout << "FIND ELEMENT" << endl;
        int size = multiset.count(8263671LL);
        chrono::high_resolution_clock::time_point t5 = chrono::high_resolution_clock::now();
        auto tuples = multiset.equal_range(8263671LL);
        for(auto it = tuples.first; it != tuples.second; ++it){
            cout << get<0>(it->second) << "\t" << get<1>(it->second) << "\t" << get<2>(it->second) << endl;
        }
        chrono::high_resolution_clock::time_point t6 = chrono::high_resolution_clock::now();

        /*
        //vector<index> ind = indexTable(v, 5, 15);

        //uint64_t  maks3 = 0;
        //for (auto &i : ind) {
        //    maks3 = maks3 < get<0>(i) ? get<0>(i) : maks3;
        //    std::cout <<setw(10) << std::get<0>(i) << " " << setw(offset_width)<< std::get<0>(std::get<1>(i)) << " " << setw(offset_width) << std::get<1>(std::get<1>(i)) << " " << std::get<2>(std::get<1>(i)) << endl;
        //}


        auto duration1 = chrono::duration_cast<chrono::microseconds>( t2 - t1 ).count();
        auto duration2 = chrono::duration_cast<chrono::microseconds>( t4 - t3 ).count();
        auto duration3 = chrono::duration_cast<chrono::microseconds>( t6 - t5 ).count();

        cout << "DUR1 " << duration1 << endl;
        cout << "DUR2 " << duration2 << endl;
        cout << "DUR3 " << duration3 << endl;

        cout << "MAKSIMUMI" << endl;
        cout << maks1 << endl;
        cout << maks2 << endl;
        //cout << maks3 << endl;



       vector<string> ref;
       string s = (*(fasta_objects2[0])).get_data();
       ref.push_back(s);


        for(string s : v) {
            map_minimizers(indexSequences(ref, 5, 15), s , 5, 15, 500);
        }

        break;
    }
    */
    int counter =0;
    for (auto &fasta_object : fasta_objects) {
        minimizer *minimizer1;
        minimizer *minimizer2;
        uint32_t minimizer1_length;
        uint32_t minimizer2_length;
        string s = fasta_object->get_data();
        uint32_t w = 5;
        uint32_t  k = 15;
        find_minimizers3(s.c_str(), (uint16_t) s.length(), w, k, &minimizer1, s.length() - k - w + 2,
                         &minimizer1_length);
        find_minimizers4(s.c_str(), (uint16_t) s.length(), w, k, &minimizer2, s.length() - k - w + 2,
                         &minimizer2_length);
        printf("%u, %u\n", minimizer1_length, minimizer2_length);
        for (int i = 0; i < minimizer1_length; i++) {
            minimizer tmp1 = minimizer1[i];
            minimizer tmp2 = minimizer2[i];
            bool condition = tmp1.hash == tmp2.hash && tmp1.index == tmp2.index && tmp1.rev == tmp2.rev;
            //printf("==> %s\n", condition ? "TRUE" : "FALSE");
            if(!condition){
                printf("%d -> %d\n", counter, i);
                break;
            }
        }
        counter++;
    }

    return 0;
}

