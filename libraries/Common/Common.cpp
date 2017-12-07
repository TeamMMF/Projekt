#include <string>
#include <vector>
#include <tuple>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <functional>
#include <iomanip>
#include <cmath>
#include <climits>

#include "Common.hpp"


using namespace std;

char complement(char c) {
    switch (c) {
        case 'A':
            return 'T';
        case 'T':
            return 'A';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        default:
            throw runtime_error("Invalid value.");
    }
}

string find_reverse_complement(string s){
    const unsigned long length = s.length();

    string reversed(s);
    reverse(reversed.begin(), reversed.end());

    const char *nucleotides = reversed.c_str();
    char *reverse_complement = new char[length + 1];        //MEM

    for(int i = 0; i < length; i++){
        reverse_complement[i] = complement(nucleotides[i]);
    }
    reverse_complement[length] = '\0';

    string ret_val = string(reverse_complement);
    return ret_val;
}

vector<string> find_kmer(int k, string s) {
    vector<string> v;
    v.reserve(s.length() - k + 1);
    for (int i = 0, len = s.length() - k; i <= len; i++) {
        v.push_back(s.substr(i, k));
    }
    return v;
}

uint64_t minimizer_hash(string s) {

    uint64_t hash = 0;
    const char* nucleotides = s.c_str();
    const unsigned long k = s.length();

    for(int i = 0; i < k; i++){
        hash += find_hash_value(nucleotides[i]) * pow(4, k - i - 1);
    }


    return hash;
}

//Thomas Wang's integer hash function
uint64_t invertible_minimizer_hash(uint64_t x){
    uint64_t  key = x;

    key = (~key) + (key << 21); // key = (key << 21) - key - 1;
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8); // key * 265
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4); // key * 21
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key;
}

uint64_t invertible_minimizer_hash_inverse(uint64_t key){
    uint64_t tmp;

    // Invert key = key + (key << 31)
    tmp = key-(key<<31);
    key = key-(tmp<<31);

    // Invert key = key ^ (key >> 28)
    tmp = key^key>>28;
    key = key^tmp>>28;

    // Invert key *= 21
    key *= 14933078535860113213u;

    // Invert key = key ^ (key >> 14)
    tmp = key^key>>14;
    tmp = key^tmp>>14;
    tmp = key^tmp>>14;
    key = key^tmp>>14;

    // Invert key *= 265
    key *= 15244667743933553977u;

    // Invert key = key ^ (key >> 24)
    tmp = key^key>>24;
    key = key^tmp>>24;

    // Invert key = (~key) + (key << 21)
    tmp = ~key;
    tmp = ~(key-(tmp<<21));
    tmp = ~(key-(tmp<<21));
    key = ~(key-(tmp<<21));

    return key;
}

size_t find_hash_value(char c){
    switch (c) {
        case 'A':
            return 0;

        case 'C':
            return 1;

        case 'G':
            return 2;

        case 'T':
            return 3;

        default:
            throw runtime_error("Invalid value");
    }
}


std::vector<triplet> find_minimizers(int w, int k, string s) {

    vector<string> kmers = find_kmer(k, s);
    vector<triplet> minimizers;
    const unsigned long size = kmers.size();
    minimizers.reserve(size - w + 1);

    auto *hashed_kmers = new pair<uint64_t, string>[size];

    for (int i = 0; i < size; i++) {
        hashed_kmers[i] = make_pair(minimizer_hash(kmers[i]), kmers[i]);
    }

    for (int i = 0, len = size - w + 1; i < len; i++) {      //prvi prozor cijeli hashirat, sljedeci samo zadnju ktorku
        auto min = hashed_kmers[i];
        int minJ = i;
        for (int j = 1; j < w; j++) {
            auto current = hashed_kmers[i + j];
            if (min.first <= current.first) {
                continue;
            }
            min = current;
            minJ = i + j;
        }
        minimizers.emplace_back(make_tuple("asfda", minJ, min.second));
    }

    delete[] hashed_kmers;

    return minimizers;
}


std::vector<tuple<uint64_t, int, int>> find_minimizers2(int w, int k, string s){

    vector<string> kmers = find_kmer(k, s);

    const unsigned long size = kmers.size();
    const unsigned long outer_loop_length = s.length() - w - k + 1;

    vector<tuple<uint64_t, int, int>> minimizers;
    minimizers.reserve(outer_loop_length + 1);
    uint64_t* hash_buffer = new uint64_t[size];
    uint64_t* r_hash_buffer = new uint64_t[size];

    for(int i = 0; i < w; i++){
        hash_buffer[i] = invertible_minimizer_hash(minimizer_hash(kmers[i]));                                //HASH
        r_hash_buffer[i] = invertible_minimizer_hash(minimizer_hash(find_reverse_complement(kmers[i])));     //HASH
    }

    for(int i = 0; i <= outer_loop_length; i++){
        uint64_t m = UINT64_MAX;

        for(int j = 0; j < w ; j++){
            uint64_t u = hash_buffer[i+j];
            uint64_t v = r_hash_buffer[i+j];
            if(u == v){
                continue;
            }
            m = min(m, min(u,v));
        }

        for(int j = 0; j < w; j++){
            uint64_t u = hash_buffer[i+j];
            uint64_t v = r_hash_buffer[i+j];

            if(u < v && u == m) {
                if(minimizers.empty()){
                    minimizers.push_back(make_tuple(m, i+j, 0));
                }

                else {
                    tuple<uint64_t, int, int> last = minimizers.back();
                    if (get<0>(last) != m && get<1>(last) != i + j) {
                        minimizers.push_back(make_tuple(m, i + j, 0));
                    }
                }

            }
            else if(v < u && v == m){
                if(minimizers.empty()){
                    minimizers.push_back(make_tuple(m, i+j, 1));
                }

                else {
                    tuple<uint64_t, int, int> last = minimizers.back();
                    if (get<0>(last) != m && get<1>(last) != i + j) {
                        minimizers.push_back(make_tuple(m, i + j, 1));
                    }
                }
            }
        }
        int next_end = i + w;
        if(next_end < size){
            hash_buffer[next_end] = invertible_minimizer_hash(minimizer_hash(kmers[next_end]));          //HASH
            r_hash_buffer[next_end] = invertible_minimizer_hash(minimizer_hash(find_reverse_complement(kmers[next_end])));   //HASH
        }
    }

    delete[] hash_buffer;
    delete[] r_hash_buffer;

    return minimizers;
}

size_t no_hash(uint64_t x){
    return  x;
}


std::unordered_multimap<uint64_t, tuple<string, int, int>, function<size_t(uint64_t)>> indexSequence(vector<string> sequences, int w, int k){
    unordered_multimap<uint64_t, tuple<string, int, int>, function<size_t(uint64_t)>> indexTable(1000, no_hash);

    int counter = 1;
    for(string seq : sequences){
        vector<tuple<uint64_t, int, int>> minimizers = find_minimizers2(w, k , seq);
        for(tuple<uint64_t, int, int> m : minimizers){
            indexTable.emplace(get<0>(m), make_tuple("seq" + std::to_string(counter), get<1>(m), get<2>(m))); //TESTIRATI S insertom!
        }
        counter++;
    }

    return indexTable;
}

bool hit_comparator(const minimizer_hit a,
                    const minimizer_hit b){
    if(get<0>(a) < get<0>(b)) return true;
    if(get<0>(a) > get<0>(b)) return false;

    if(get<1>(a) < get<1>(b)) return true;
    if(get<1>(a) > get<1>(b)) return false;

    if(get<2>(a) < get<2>(b)) return true;
    if(get<2>(a) > get<2>(b)) return false;

    if(get<3>(a) < get<3>(b)) return true;
    if(get<3>(a) > get<3>(b)) return false;

    return true;
}


std::vector<hashMinPair> indexTable(vector<string> sequences, int w, int k){
    vector<hashMinPair> table;

    int counter = 1;
    for(string seq : sequences){
        vector<tuple<uint64_t, int, int>> minimizers = find_minimizers2(w, k , seq);
        for(tuple<uint64_t, int, int> m : minimizers){
            table.emplace_back(get<0>(m), make_tuple("seq" + std::to_string(counter), get<1>(m), get<2>(m))); //TESTIRATI S insertom!
        }
        counter++;
    }

    std::sort(table.begin(), table.end());
    return table;
}




vector<mapInfo> map_minimizers(unordered_multimap<uint64_t, tuple<string, int, int>, function<size_t(uint64_t)>> lookup_table,
         string query_sequence,
         int w,
         int k,
         int epsilon){


    vector<minimizer_hit> hits;
    auto minimizers = find_minimizers2(w, k, query_sequence);


    for(auto minimizer : minimizers) {
        uint64_t h = get<0>(minimizer);
        int i_q = get<1>(minimizer);
        int r_q = get<2>(minimizer);
        auto tirs = lookup_table.equal_range(h);

        for (auto it = tirs.first; it != tirs.second; ++it) {
            auto tir = it->second;
            string t = get<0>(tir);
            int i2 = get<1>(tir);
            int r2 = get<2>(tir);
            int same_strand = r_q == r2 ? 0 : 1;

            if(r_q == 0) {
                hits.emplace_back(make_tuple(t, 0, i_q - i2, i2));
            }
            else {
                hits.emplace_back(make_tuple(t, 1, i_q + i2, i2));
            }

        }
    }

    sort(hits.begin(),hits.end(),hit_comparator);

    vector<mapInfo> info;

    int b = 1;
    for(int e = 1, limit = hits.size(); e<limit; e++) {
        string t1 = get<0>(hits[e - 1]);
        string t2 = get<0>(hits[e]);
        int r1 = get<1>(hits[e - 1]);
        int r2 = get<1>(hits[e]);
        int diff1 = get<2>(hits[e - 1]);
        int diff2 = get<2>(hits[e]);
        if (r1 != r2 || t1 != t2 || diff2 - diff1 >= epsilon) {

            int min_i_prime = INT_MAX;
            int max_i_prime = 0;
            int min_i = INT_MAX;
            int max_i = 0;
            for(int i = b; i <= e; i++){
                int i_prime_tmp = get<3>(hits[i]);
                int i_tmp = get<1>(hits[i]) ? get<2>(hits[i]) - i_prime_tmp : get<2>(hits[i]) + i_prime_tmp;
                if(i_prime_tmp < min_i_prime){
                    min_i_prime = i_prime_tmp;
                    min_i = i_tmp;
                }

                if(i_prime_tmp > max_i_prime){
                    max_i_prime = i_prime_tmp;
                    max_i = i_tmp;
                }

            }
            if(max_i - min_i <100){
                continue;
            }
            mapInfo mi;
            mi.query_min_index = min_i;
            mi.query_max_index = max_i;
            mi.reverse = r1 != 0 ? true : false;
            mi.target_min_index = min_i_prime;
            mi.target_max_index = max_i_prime;

            info.emplace_back(mi);

            b = e + 1;

            /*
            cout << mi.target_min_index << " , "
                 << mi.target_max_index << " , "
                 << mi.reverse << " , "
                 << mi.query_min_index << " , "
                 << mi.query_max_index << endl;
            */
        }

    }

    return info;
}

