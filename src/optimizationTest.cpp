//
// Created by matep on 08/12/2017.
//

#include <vector>
#include <tuple>
#include <chrono>
#include <iostream>
#include <algorithm>
#include <CustomTypes.h>
#include "Common.hpp"

using namespace std;

int main(){

    std::string s = "AACCTTGGAACCGGTTACGTGCTAGCAGTGATGCTGAGCTGAGAGATCTTAGAGCTAGTCAGCTACGATCAGCTACGCTACGACTACGATTATTAAGCGGGCGGGATCACGACTACGACTAGCGACTTATGGAGTCTCTCTTATTAGGGTGTGGTTCTCTGCGCGTATAGGCTGATCGATCAGCTAGGTGAGCTAGCATCGATCAGTG";
    /*
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


//******************************************************************************************

    // NE ZNAM ZAŠTO RADIŠ OVO, BUT PLS NE MIJENJAJ FORMATE ISPISA VIŠE


//******************************************************************************************


    printf("%lu\n", minimizer_hash(s));
    printf("%lu\n", minimizer_hash3(s.c_str(), s.length()));

    uint16_t w = 5;
    minimizer *minimizers;
    uint32_t actual_length;
    find_minimizers3(s.c_str(),(uint16_t) s.length(), w, k, &minimizers, s.length() - k - w + 2, &actual_length);
    //for(int i = 0; i < actual_length; i++){
    //    printf("(%lu, %d, %s)\n", minimizers[i].hash, minimizers[i].index, minimizers[i].rev ? "true" : "false");
    //}
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
    //for(auto m : minis){
    //    printf("(%ld, %d, %s)\n", get<0>(m), get<1>(m), get<2>(m) != 0 ? "true" : "false");
    //}

    uint64_t stl_time = 0;
    uint64_t c_time = 0;
    uint64_t c_red_time = 0;


    for(int i = 0; i < 100000; i++){
        minimizer* minimizer1;
        minimizer* minimizer2;
        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
        find_minimizers2(w,k,s);
        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        find_minimizers3(s.c_str(),(uint16_t) s.length(), w, k, &minimizer1, s.length() - k - w + 2, &actual_length);
        chrono::high_resolution_clock::time_point t3 = chrono::high_resolution_clock::now();
        find_minimizers4(s.c_str(),(uint16_t) s.length(), w, k, &minimizer2, s.length() - k - w + 2, &actual_length);
        chrono::high_resolution_clock::time_point t4 = chrono::high_resolution_clock::now();
        auto duration1 = chrono::duration_cast<chrono::microseconds>( t2 - t1 ).count();
        auto duration2 = chrono::duration_cast<chrono::microseconds>( t3 - t2 ).count();
        auto duration3 = chrono::duration_cast<chrono::microseconds>( t4 - t3 ).count();
        stl_time += duration1;
        c_time += duration2;
        c_red_time += duration3;
        free(minimizer1);
        free(minimizer2);
    }

    //stl_time /= 1000;
    //c_time /= 1000;
    printf("STL  time: %lu\n", stl_time);
    printf("C    time: %lu\n", c_time);
    printf("Cred time: %lu\n", c_red_time);

    /*
    std::unordered_multimap<uint64_t, int> hash_to_index_map_addr;
    minimizer* ordered_minimizers_addr;
    uint32_t ordered_minimizers_addr_l;

    process_sequence(s.c_str(), s.length(), w, k, hash_to_index_map_addr, &ordered_minimizers_addr, &ordered_minimizers_addr_l);

    printf("PROCESS SEQUENCE: \n");
    for(int i = 0; i < ordered_minimizers_addr_l; i++){
        minimizer tmp = ordered_minimizers_addr[i];
        printf("(%lu, %d, %s)\n", tmp.hash, tmp.index, tmp.rev ? "True" : "False");
    }
    printf("%u\n", ordered_minimizers_addr_l);
    printf("MAP INFO\n");
    std::pair <std::unordered_multimap<uint64_t ,int>::iterator, std::unordered_multimap<uint64_t ,int>::iterator> range;
    for(auto i = hash_to_index_map_addr.begin(); i != hash_to_index_map_addr.end(); i = range.second)
    {
        // Get the range of the current key
        range = hash_to_index_map_addr.equal_range(i->first);

        // Now print out that whole range
        for(auto d = range.first; d != range.second; ++d)
            std::cout << d->first << ": " << d->second << '\n';
    }

    printf("%ld", hash_to_index_map_addr.size());
     */
    /*
    free(minimizers);
    minimizer* minimizer1;
    minimizer* minimizer2;
    uint32_t minimizer1_length;
    uint32_t minimizer2_length;

    find_minimizers3(s.c_str(),(uint16_t) s.length(), w, k, &minimizer1, s.length() - k - w + 2, &minimizer1_length);
    find_minimizers4(s.c_str(),(uint16_t) s.length(), w, k, &minimizer2, s.length() - k - w + 2, &minimizer2_length);
    printf("%u, %u\n" ,minimizer1_length, minimizer2_length);
    for(int i = 0; i < minimizer1_length; i++){
        minimizer tmp1 = minimizer1[i];
        minimizer tmp2 = minimizer2[i];
        bool condition = tmp1.hash == tmp2.hash && tmp1.index == tmp2.index && tmp1.rev == tmp2.rev;
        printf("(%lu, %d, %s), (%lu, %d, %s) ==> %s\n", tmp1.hash, tmp1.index, tmp1.rev ? "True" : "False", tmp2.hash, tmp2.index, tmp2.rev ? "True" : "False", condition ? "TRUE" : "FALSE");

    }

    free(minimizer1);
    free(minimizer2);
    */

    printf("%llu\n", invertible_minimizer_hash_inverse(0));


    string s1 = "GAGATCTGA";
    string s2 = "CAGAACTGT";
    string s3 = "GAGTTCAGA";
    /*
    std::vector<vector<uint64_t>> minimizers(3);
    std::vector<hashMinPair> minimizer_pairs;

    find_minimizers5(s1.c_str(), s1.length(), 0, 2 ,3, &(minimizers[0]), &minimizer_pairs);
    find_minimizers5(s2.c_str(), s2.length(), 1, 2 ,3, &(minimizers[1]), &minimizer_pairs);
    find_minimizers5(s3.c_str(), s3.length(), 2, 2 ,3, &(minimizers[2]), &minimizer_pairs);


        for(int j = 0; j < minimizer_pairs.size(); j++){
            printf("%llu, %u, %u, %s\n", minimizer_pairs[j].hash, minimizer_pairs[j].seq_id, minimizer_pairs[j].index, minimizer_pairs[j].rev ? "True" : "False");
        }

    sort(minimizer_pairs.begin(), minimizer_pairs.end(), hashMinPair_comparator);
    printf("\n");
    for(int j = 0; j < minimizer_pairs.size(); j++){
        printf("%llu, %u, %u, %s\n", minimizer_pairs[j].hash, minimizer_pairs[j].seq_id, minimizer_pairs[j].index, minimizer_pairs[j].rev ? "True" : "False");
    }

    printf("%d", sizeof(std::tuple<uint64_t,uint32_t ,uint32_t, bool>));
     */

    std::vector<minimizer> first;
    std::vector<minimizer> second;
    std::vector<minimizer> third;
    std::unordered_map<uint64_t, uint32_t> map1;
    std::unordered_map<uint64_t, hashMinPair2> map2;
    find_minimizers7(s.c_str(), s.length(), 0, 3, 4, first);
    find_minimizers_deq(s.c_str(), s.length(), 0, 3, 4, second, map1, map2  );
    find_minimizers_deq_single(s.c_str(), s.length(), 0, 3, 4, third);

    printf("%lu -> %lu -> %lu\n", first.size(), second.size(), third.size());
    for(int i = 0; i < third.size(); i++){
        printf("(%22lu, %10d, %5s) -> (%22lu, %10d, %5s) -> (%22lu, %10d, %5s)\n", first[i].hash, first[i].index, first[i].rev ? "True" : "False", second[i].hash, second[i].index, second[i].rev ? "True" : "False", third[i].hash, third[i].index, third[i].rev ? "True" : "False" );
    }


    uint64_t old_time = 0;
    uint64_t ddeq_time = 0;
    uint64_t sdeq_time = 0;
    for(int i = 0; i < 10000; i++){
        std::vector<minimizer> first;
        std::vector<minimizer> second;
        std::vector<minimizer> third;
        chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
        find_minimizers7(s.c_str(), s.length(), 0, 5, 15, first);
        chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
        find_minimizers_deq(s.c_str(), s.length(), 0, 5, 15, second, map1, map2  );
        chrono::high_resolution_clock::time_point t3 = chrono::high_resolution_clock::now();
        find_minimizers_deq_single(s.c_str(), s.length(), 0, 5, 15, third);
        chrono::high_resolution_clock::time_point t4 = chrono::high_resolution_clock::now();

        auto duration1 = chrono::duration_cast<chrono::microseconds>( t2 - t1 ).count();
        auto duration2 = chrono::duration_cast<chrono::microseconds>( t3 - t2 ).count();
        auto duration3 = chrono::duration_cast<chrono::microseconds>( t4 - t3 ).count();
        old_time += duration1;
        ddeq_time += duration2;
        sdeq_time += duration3;
    }

    printf("Staro vrijeme: %22lu\n", old_time / 10000);
    printf("DDEQ vrijeme:  %22lu\n", ddeq_time / 10000);
    printf("SDEQ vrijeme:  %22lu\n", sdeq_time / 10000);

    return 0;
}


