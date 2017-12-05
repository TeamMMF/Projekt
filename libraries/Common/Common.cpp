#include <string>
#include <vector>
#include <tuple>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <functional>
#include <iomanip>
#include <math.h>
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
    delete(reverse_complement);
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

    delete(nucleotides);

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

    delete(hashed_kmers);

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
        hash_buffer[i] = minimizer_hash(kmers[i]);
        r_hash_buffer[i] = minimizer_hash(find_reverse_complement(kmers[i]));
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
            hash_buffer[next_end] = minimizer_hash(kmers[next_end]);
            r_hash_buffer[next_end] = minimizer_hash(find_reverse_complement(kmers[next_end]));
        }
    }

    delete(hash_buffer);
    delete(r_hash_buffer);

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

