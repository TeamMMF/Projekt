//
// Created by matep on 08/12/2017.
//

#include <vector>
#include <tuple>
#include <chrono>
#include "Common.hpp"

using namespace std;

int main(){
    std::string s = "AACCTTGGAACCGGTTACGTGCTAGCAGTGATGCTGAGCTGAGAGATCTTAGAGCTAGTCAGCTACGATCAGCTACGCTACGACTACGATTATTAAGCGGGCGGGATCACGACTACGACTAGCGACTTATGGAGTCTCTCTTATTAGGGTGTGGTTCTCTGCGCGTATAGGCTGATCGATCAGCTAGGTGAGCTAGCATCGATCAGTG";
    uint8_t k = 3;
    char** kmers;
    find_kmers(s.c_str(), 3, &kmers, s.length() - k + 1);
    printf("VRATIO SE IZ FUNKCIJE\n");
    printf("VRATIO SE\n");
    //for(int i = 0; i < s.length() - k + 1 ; i++){
    //
    // }

    destroy_kmers(&kmers, (uint16_t) (s.length() - k + 1));

    char* rev_comp;
    find_reverse_complement(s.c_str(), s.length(), &rev_comp);
    printf("%s\n", s.c_str());
    printf("%s\n", rev_comp);

    destroy_reverse_complement(&rev_comp);

    printf("%ld\n", minimizer_hash(s));
    printf("%ld\n", minimizer_hash3(s.c_str(), s.length()));

    uint16_t w = 5;
    minimizer *minimizers;
    uint32_t actual_length;
    find_minimizers3(s.c_str(),(uint16_t) s.length(), w, k, &minimizers, s.length() - k - w + 2, &actual_length);
    for(int i = 0; i < actual_length; i++){
        printf("(%ld, %d, %s)\n", minimizers[i].hash, minimizers[i].index, minimizers[i].rev ? "true" : "false");
    }
    printf("\n\n");
    vector<tuple<uint64_t, int, int>> minis = find_minimizers2(w, k, s);
    printf("\n\n");
    char k10[] = {'C', 'C', 'G'};
    char k11[] = {'C', 'G', 'G'};
    char k12[] = {'G', 'G', 'T'};
    char k13[] = {'G', 'T', 'T'};
    printf("%22lu\n", invertible_minimizer_hash(minimizer_hash3(k10, k)));
    printf("%22lu\n", invertible_minimizer_hash(minimizer_hash3(k11, k)));
    printf("%22lu\n", invertible_minimizer_hash(minimizer_hash3(k12, k)));
    printf("%22lu\n", invertible_minimizer_hash(minimizer_hash3(k13, k)));
    printf("\n");
    for(auto m : minis){
        printf("(%ld, %d, %s)\n", get<0>(m), get<1>(m), get<2>(m) != 0 ? "true" : "false");
    }

    uint64_t stl_time = 0;
    uint64_t c_time = 0;

    for(int i = 0; i < 100000; i++){
        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
        find_minimizers2(w,k,s);
        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        find_minimizers3(s.c_str(),(uint16_t) s.length(), w, k, &minimizers, s.length() - k - w + 2, &actual_length);
        chrono::high_resolution_clock::time_point t3 = chrono::high_resolution_clock::now();
        auto duration1 = chrono::duration_cast<chrono::microseconds>( t2 - t1 ).count();
        auto duration2 = chrono::duration_cast<chrono::microseconds>( t3 - t2 ).count();
        stl_time += duration1;
        c_time += duration2;
    }
    //stl_time /= 1000;
    //c_time /= 1000;
    printf("STL time: %lu\n", stl_time);
    printf("C   time: %lu\n", c_time);
}
