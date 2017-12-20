

#include "Common.hpp"
#include "../threadpool/include/thread_pool/thread_pool.hpp"
#include <tuple>
#include <algorithm>
#include <unordered_set>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <climits>
#include <cstring>
#include <CustomTypes.h>


using namespace std;

string find_reverse_complement(string s) {
    const unsigned long length = s.length();

    string reversed(s);
    reverse(reversed.begin(), reversed.end());

    const char *nucleotides = reversed.c_str();
    char *reverse_complement = new char[length + 1];        //MEM

    for (int i = 0; i < length; i++) {
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
    const char *nucleotides = s.c_str();
    const unsigned long k = s.length();

    for (int i = 0; i < k; i++) {
        hash += find_hash_value(nucleotides[i]) * pow(4, k - i - 1);
    }


    return hash;
}

//Thomas Wang's integer hash function
uint64_t invertible_minimizer_hash(uint64_t x) {
    uint64_t key = x;

    key = (~key) + (key << 21); // key = (key << 21) - key - 1;
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8); // key * 265
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4); // key * 21
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key;
}

uint64_t invertible_minimizer_hash_inverse(uint64_t key) {
    uint64_t tmp;

    // Invert key = key + (key << 31)
    tmp = key - (key << 31);
    key = key - (tmp << 31);

    // Invert key = key ^ (key >> 28)
    tmp = key ^ key >> 28;
    key = key ^ tmp >> 28;

    // Invert key *= 21
    key *= 14933078535860113213u;

    // Invert key = key ^ (key >> 14)
    tmp = key ^ key >> 14;
    tmp = key ^ tmp >> 14;
    tmp = key ^ tmp >> 14;
    key = key ^ tmp >> 14;

    // Invert key *= 265
    key *= 15244667743933553977u;

    // Invert key = key ^ (key >> 24)
    tmp = key ^ key >> 24;
    key = key ^ tmp >> 24;

    // Invert key = (~key) + (key << 21)
    tmp = ~key;
    tmp = ~(key - (tmp << 21));
    tmp = ~(key - (tmp << 21));
    key = ~(key - (tmp << 21));

    return key;
}

size_t find_hash_value(char c) {
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


std::vector<tuple<uint64_t, int, int>> find_minimizers2(int w, int k, string s) {

    vector<string> kmers = find_kmer(k, s);

    const unsigned long size = kmers.size();
    const unsigned long outer_loop_length = s.length() - w - k + 1;

    vector<tuple<uint64_t, int, int>> minimizers;
    minimizers.reserve(outer_loop_length + 1);
    uint64_t *hash_buffer = new uint64_t[size];
    uint64_t *r_hash_buffer = new uint64_t[size];

    for (int i = 0; i < w; i++) {
        hash_buffer[i] = invertible_minimizer_hash(minimizer_hash(kmers[i]));                                //HASH
        r_hash_buffer[i] = invertible_minimizer_hash(minimizer_hash(find_reverse_complement(kmers[i])));     //HASH
    }

    for (int i = 0; i <= outer_loop_length; i++) {
        uint64_t m = UINT64_MAX;

        for (int j = 0; j < w; j++) {
            uint64_t u = hash_buffer[i + j];
            uint64_t v = r_hash_buffer[i + j];
            if (u == v) {
                continue;
            }
            m = min(m, min(u, v));
        }

        for (int j = 0; j < w; j++) {
            uint64_t u = hash_buffer[i + j];
            uint64_t v = r_hash_buffer[i + j];

            if (u < v && u == m) {
                if (minimizers.empty()) {
                    minimizers.push_back(make_tuple(m, i + j, 0));
                } else {
                    tuple<uint64_t, int, int> last = minimizers.back();
                    if (get<0>(last) != m && get<1>(last) != i + j) {
                        minimizers.push_back(make_tuple(m, i + j, 0));
                    }
                }

            } else if (v < u && v == m) {
                if (minimizers.empty()) {
                    minimizers.push_back(make_tuple(m, i + j, 1));
                } else {
                    tuple<uint64_t, int, int> last = minimizers.back();
                    if (get<0>(last) != m && get<1>(last) != i + j) {
                        minimizers.push_back(make_tuple(m, i + j, 1));
                    }
                }
            }
        }
        int next_end = i + w;
        if (next_end < size) {
            hash_buffer[next_end] = invertible_minimizer_hash(minimizer_hash(kmers[next_end]));          //HASH
            r_hash_buffer[next_end] = invertible_minimizer_hash(
                    minimizer_hash(find_reverse_complement(kmers[next_end])));   //HASH
        }
    }

//    for(int i = 0; i < size; i++){
//        printf("%2d h=%22llu, hr = %22llu\n", i, hash_buffer[i], r_hash_buffer[i]);
//    }

    delete[] hash_buffer;
    delete[] r_hash_buffer;

    return minimizers;
}


size_t no_hash(uint64_t x) {
    return x;
}


std::unordered_multimap<uint64_t, tuple<string, int, int>, function<size_t(uint64_t)>> indexSequences(vector<string> sequences, int w, int k) {

    unordered_multimap<uint64_t, tuple<string, int, int>, function<size_t(uint64_t)>> indexTable(1000, no_hash);

    int counter = 1;
    for (string seq : sequences) {
        vector<tuple<uint64_t, int, int>> minimizers = find_minimizers2(w, k, seq);
        for (tuple<uint64_t, int, int> m : minimizers) {
            indexTable.emplace(get<0>(m), make_tuple("seq" + std::to_string(counter), get<1>(m),
                                                     get<2>(m))); //TESTIRATI S insertom!
        }
        counter++;
    }

    return indexTable;
}

std::unordered_multimap<uint64_t, hashEntry, function<size_t(uint64_t)>> indexSequence(string sequence, int w, int k) {

    //no hash??
    unordered_multimap<uint64_t, hashEntry, function<size_t(uint64_t)>> indexTable(1000, no_hash);

    minimizer** minimizers;
    uint32_t  min_l_real, min_l_pred;

    //sto je min_l_pred?
    find_minimizers3(sequence.c_str(), sequence.length(), w, k, minimizers, min_l_pred, &min_l_real);

    for (int i = 0; i < min_l_real; i++) {
        minimizer *m = minimizers[i];

        hashEntry *entry;
        entry->index = m->index;
        entry->rev = m-> rev;

        //todo napraviti ovo bez errora
        // indexTable.emplace(m->hash, entry);
    }

    return  indexTable;
}

bool hit_comparator(const minimizer_hit a,
                    const minimizer_hit b) {
    if (get<0>(a) < get<0>(b)) return true;
    if (get<0>(a) > get<0>(b)) return false;

    if (get<1>(a) < get<1>(b)) return true;
    if (get<1>(a) > get<1>(b)) return false;

    if (get<2>(a) < get<2>(b)) return true;
    if (get<2>(a) > get<2>(b)) return false;

    if (get<3>(a) < get<3>(b)) return true;
    if (get<3>(a) > get<3>(b)) return false;

    return true;
}


std::vector<hashMinPair3> indexTable(vector<string> sequences, int w, int k) {
    vector<hashMinPair3> table;

    int counter = 1;
    for (string seq : sequences) {
        vector<tuple<uint64_t, int, int>> minimizers = find_minimizers2(w, k, seq);
        for (tuple<uint64_t, int, int> m : minimizers) {
            table.emplace_back(get<0>(m), make_tuple("seq" + std::to_string(counter), get<1>(m),
                                                     get<2>(m))); //TESTIRATI S insertom!
        }
        counter++;
    }

    std::sort(table.begin(), table.end());
    return table;
}


vector<mapInfo>  map_minimizers(unordered_multimap<uint64_t, tuple<string, int, int>, function<size_t(uint64_t)>> lookup_table,
                                string query_sequence,
                                int w,
                                int k,
                                int epsilon) {


    vector<minimizer_hit> hits;
    auto minimizers = find_minimizers2(w, k, query_sequence);


    for (auto minimizer : minimizers) {
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

            if (r_q == 0) {
                hits.emplace_back(make_tuple(t, 0, i_q - i2, i2));
            } else {
                hits.emplace_back(make_tuple(t, 1, i_q + i2, i2));
            }

        }
    }

    sort(hits.begin(), hits.end(), hit_comparator);

    vector<mapInfo> info;

    int b = 1;
    for (int e = 1, limit = hits.size(); e < limit; e++) {
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
            for (int i = b; i <= e; i++) {
                int i_prime_tmp = get<3>(hits[i]);
                int i_tmp = get<1>(hits[i]) ? get<2>(hits[i]) - i_prime_tmp : get<2>(hits[i]) + i_prime_tmp;
                if (i_prime_tmp < min_i_prime) {
                    min_i_prime = i_prime_tmp;
                    min_i = i_tmp;
                }

                if (i_prime_tmp > max_i_prime) {
                    max_i_prime = i_prime_tmp;
                    max_i = i_tmp;
                }

            }
            if (max_i - min_i < 100) {
                continue;
            }
            mapInfo mi;
            mi.query_min_index = min_i;
            mi.query_max_index = max_i;
            mi.reverse = r1 != 0;
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

/************************************************* NON STL C++ CODE ********************************************************/

/**
 * Returns the complement base
 * @param c
 * @return
 */
char complement(char c)
{
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

/**
 * Calculates the reverse complement od the given char array.
 * CALL FREE ON RETURNED POINTER
 * @param seq
 * @param seq_l - length without null terminator
 * @param rev_comp
 * @return
 */
void find_reverse_complement(const char *seq, uint32_t seq_l, char** rev_comp)
{
    *rev_comp = (char*) malloc((seq_l + 1)* sizeof(char));

    uint16_t counter = 0;
    for (int i = seq_l - 1; i >= 0; i--) {
        (*rev_comp)[counter++] = complement(seq[i]);
    }
    (*rev_comp)[seq_l] = '\0';

    return;
}
/**
 * Frees the malloc for the reverse_complement
 * @param rev_comp
 */
void destroy_reverse_complement(char** rev_comp)
{
    free(*rev_comp);
}

/**
 * Finds all kmers of the given sequnce
 * CALL DESTROY_KMERS
 * CHECK BEFORE CALLING THATH seq_l -k + 1 > 0
 * @param seq
 * @param k
 * @param kmers
 * @param kmers_l
 */
void find_kmers(const char *seq, uint32_t k, char ***kmers, uint32_t kmers_l)
{
    *kmers = (char**) malloc(kmers_l*sizeof(char*));
    for(int i = 0; i < kmers_l; i++){
        (*kmers)[i] = (char *) malloc((k+1)*sizeof(char));
        if((*kmers)[i] == nullptr){
            exit(1);
        }
        strncpy((*kmers)[i], &(seq[i]), k);
        (*kmers)[i][k] = '\0';
        //printf("%s %p\n", (*kmers)[i], &(*kmers)[i]);
    }

    return;
}

/**
 * Frees all the memory that was given to the given kmer array
 * @param kmers
 * @param kmers_l
 */
void destroy_kmers(char*** kmers, uint32_t kmers_l)
{
    for(int i = 0; i < kmers_l; i++){
        //printf("Oslobadam: %p\n", &(*kmers)[i]);
        free((*kmers)[i]);
    }
    free(*kmers);
}

/**
 * Calculates hash value of given char array
 * @param seq
 * @param seq_l
 * @return
 */
uint64_t minimizer_hash3(const char* seq, uint32_t seq_l)
{
    uint64_t hash = 0;

    for (int i = 0; i < seq_l; i++) {
        hash += find_hash_value(seq[i]) * pow(4, seq_l - i - 1);
    }

    return hash;
}

uint64_t minimizer_hash3_rev(const char* seq, uint32_t seq_l)
{
    uint64_t hash = 0;

    for (int i = seq_l - 1; i >= 0; i--) {
        hash += find_hash_value(complement(seq[i])) * pow(4, i);
    }

    return hash;
}


/*
 * Predati mallocirane minimizers!
 * Freeati u pocetnoj funkciji!
 */

void find_minimizers3
        (const char *seq,
         uint32_t seq_l,
         uint32_t w,
         uint32_t k,
         minimizer** minimizers,
         uint32_t min_l_pred,
         uint32_t* min_l_real
        )
{
    *minimizers = (minimizer*) malloc(min_l_pred * sizeof(minimizer));
    *min_l_real = 0;

    const uint32_t kmers_l = seq_l - k + 1;
    char **kmers = (char **) malloc(kmers_l);
    find_kmers(seq, k, &kmers, kmers_l);
    //const uint16_t min_l = seq_l - w - k + 2;

    uint64_t *hash_buffer = new uint64_t[kmers_l];
    uint64_t *r_hash_buffer = new uint64_t[kmers_l];

    for (int i = 0; i < w; i++) {
        hash_buffer[i] = invertible_minimizer_hash(minimizer_hash3(kmers[i], k));
        char* rev_comp;
        find_reverse_complement(kmers[i],k, &rev_comp);
        r_hash_buffer[i] = invertible_minimizer_hash(minimizer_hash3(rev_comp,k));   //HASH
        destroy_reverse_complement(&rev_comp);
    }

    uint16_t min_position = 0; //position where minimizer should be constructed next
    uint64_t last_min_hash = UINT64_MAX;
    int64_t last_min_position = -1;

    for (uint32_t i = 0; i < min_l_pred; i++) {
        uint64_t u;
        uint64_t v;

        if(last_min_position != -1 && last_min_position >= i && last_min_position < i + w){
            u = hash_buffer[i + w - 1];
            v = r_hash_buffer[i + w - 1];

            if(u == v){
                continue;
            }

            else if(u < v && u <= last_min_hash){
                (*minimizers)[min_position++] = (minimizer) {u, (uint32_t)(i + w - 1), false};
                last_min_position = i + w - 1;
                last_min_hash = u;
            }

            else if(u > v && v <= last_min_hash) {
                (*minimizers)[min_position++] = (minimizer) {v, (uint32_t)(i + w - 1), true};
                last_min_position = i + w - 1;
                last_min_hash = v;
            }
        }
        else {
            uint64_t m = UINT64_MAX;  //van petlje init?

            uint32_t *min_positions = new uint32_t[w];
            bool *min_rev = new bool[w];
            uint16_t min_pos_size = 0;
            for (int j = 0; j < w; j++) {
                u = hash_buffer[i + j];
                v = r_hash_buffer[i + j];

                if(u == v){
                    continue;
                }

                if(u < m || v <  m){
                    min_pos_size = 0;

                    if(u < v){
                        min_positions[min_pos_size] = i + j;
                        min_rev[min_pos_size] = false;
                        m = u;

                    } else {
                        min_positions[min_pos_size] = i + j;
                        min_rev[min_pos_size] = true;
                        m = v;
                    }

                    min_pos_size++;
                }

                else if(u == m){
                    min_positions[min_pos_size] = i + j;
                    min_rev[min_pos_size++] = false;
                }

                else if(v == m){
                    min_positions[min_pos_size] = i + j;
                    min_rev[min_pos_size++] = true;
                }
            }

            last_min_hash = m;
            last_min_position = min_positions[min_pos_size - 1];

            for(uint32_t j = 0; j < min_pos_size; j++){
                (*minimizers)[min_position++] = (minimizer) {m, (uint32_t)(min_positions[j]), min_rev[j]};
            }


            delete[] min_positions;
            delete[] min_rev;
        }

        int next_end = i + w;
        if (next_end < kmers_l) {
            hash_buffer[next_end] = invertible_minimizer_hash(minimizer_hash3(kmers[next_end],k));//HASH
            char*rev_comp;
            find_reverse_complement(kmers[next_end], k, &rev_comp);
            r_hash_buffer[next_end] = invertible_minimizer_hash(minimizer_hash3(rev_comp,k));   //HASH
            //printf("%d %s %llu, %s %llu\n",next_end, kmers[next_end], invertible_minimizer_hash(minimizer_hash3(kmers[next_end],k)),
            //                             rev_comp, invertible_minimizer_hash(minimizer_hash3(rev_comp,k)));
            destroy_reverse_complement(&rev_comp);
        }
    }

    *min_l_real = min_position;
//    for(int i = 0; i < kmers_l; i++){
//        printf("%2d h=%22llu, hr = %22llu\n", i, hash_buffer[i], r_hash_buffer[i]);
//    }
//    for(int i = 0; i < *min_l_real; i++){
//        printf("(%lld, %d, %s)\n", (*minimizers)[i].hash, (*minimizers)[i].index, (*minimizers)[i].rev ? "true" : "false");
//    }
//
//    printf("\n\n");

    //FREE BLOK
    destroy_kmers(&kmers, kmers_l);
    delete[] hash_buffer;
    delete[] r_hash_buffer;
    //

    return;
}

void find_minimizers4
        (const char *seq,
         uint32_t seq_l,
         uint32_t w,
         uint32_t k,
         minimizer** minimizers,
         uint32_t min_l_pred,
         uint32_t* min_l_real
        )
{
    uint32_t minimizers_current_size = min_l_pred / 4 + 1;
    *minimizers = (minimizer*) malloc(minimizers_current_size * sizeof(minimizer));
    *min_l_real = 0;

    const uint32_t kmers_l = seq_l - k + 1;
    //const uint16_t min_l = seq_l - w - k + 2;

    uint64_t *hash_buffer = new uint64_t[w];
    uint64_t *r_hash_buffer = new uint64_t[w];

    for (int i = 0; i < w; i++) {
        hash_buffer[i] = invertible_minimizer_hash(minimizer_hash3(&(seq[i]), k));
        r_hash_buffer[i] = invertible_minimizer_hash(minimizer_hash3_rev(&(seq[i]),k));   //HASH
    }

    uint16_t min_position = 0; //position where minimizer should be constructed next
    uint64_t last_min_hash = UINT64_MAX;
    int64_t last_min_position = -1;

    for (uint32_t i = 0; i < min_l_pred; i++) {
        uint64_t u;
        uint64_t v;

        if(last_min_position != -1 && last_min_position >= i && last_min_position < i + w){
            u = hash_buffer[(i + w - 1) % w];
            v = r_hash_buffer[(i + w - 1) % w];

            if(u == v){
                continue;
            }

            else if(u < v && u <= last_min_hash){
                if(min_position == minimizers_current_size){
                    *minimizers = (minimizer*) realloc(*minimizers, minimizers_current_size*2* sizeof(minimizer));
                    if(*minimizers == nullptr){
                        fprintf(stdout, "COULDNT REALLOC MEMORY\n");
                        return;
                    }
                    minimizers_current_size *= 2;
                    //printf("REALLOCED\n");
                }
                (*minimizers)[min_position++] = (minimizer) {u, (uint32_t)(i + w - 1), false};      //REALLOC
                last_min_position = i + w - 1;
                last_min_hash = u;
            }

            else if(u > v && v <= last_min_hash) {
                if(min_position == minimizers_current_size){
                    *minimizers = (minimizer*) realloc(*minimizers, minimizers_current_size*2* sizeof(minimizer));
                    if(*minimizers == nullptr){
                        fprintf(stdout, "COULDNT REALLOC MEMORY\n");
                        return;
                    }
                    minimizers_current_size *= 2;
                    //printf("REALLOCED\n");
                }
                (*minimizers)[min_position++] = (minimizer) {v, (uint32_t)(i + w - 1), true};       //REALLOC
                last_min_position = i + w - 1;
                last_min_hash = v;
            }
        }
        else {
            uint64_t m = UINT64_MAX;

            uint32_t *min_positions = new uint32_t[w];
            bool *min_rev = new bool[w];
            uint16_t min_pos_size = 0;
            for (int j = 0; j < w; j++) {
                u = hash_buffer[(i + j) % w];
                v = r_hash_buffer[(i + j) % w];

                if(u == v){
                    continue;
                }

                if(u < m || v <  m){
                    min_pos_size = 0;

                    if(u < v){
                        min_positions[min_pos_size] = i + j;
                        min_rev[min_pos_size] = false;
                        m = u;

                    } else {
                        min_positions[min_pos_size] = i + j;
                        min_rev[min_pos_size] = true;
                        m = v;
                    }

                    min_pos_size++;
                }

                else if(u == m){
                    min_positions[min_pos_size] = i + j;
                    min_rev[min_pos_size++] = false;
                }

                else if(v == m){
                    min_positions[min_pos_size] = i + j;
                    min_rev[min_pos_size++] = true;
                }
            }

            last_min_hash = m;
            last_min_position = min_positions[min_pos_size - 1];

            for(uint32_t j = 0; j < min_pos_size; j++){
                if(min_position == minimizers_current_size){
                    *minimizers = (minimizer*) realloc(*minimizers, minimizers_current_size*2* sizeof(minimizer));
                    if(*minimizers == nullptr){
                        fprintf(stdout, "COULDNT REALLOC MEMORY\n");
                        return;
                    }
                    minimizers_current_size *= 2;
                    //printf("REALLOCED\n");
                }
                (*minimizers)[min_position++] = (minimizer) {m, (uint32_t)(min_positions[j]), min_rev[j]}; // REALLOC
            }


            delete[] min_positions;
            delete[] min_rev;
        }

        int next_end = i + w;
        if (next_end < kmers_l) {
            hash_buffer[next_end % w] = invertible_minimizer_hash(minimizer_hash3(&(seq[next_end]),k));//HASH
            r_hash_buffer[next_end % w] = invertible_minimizer_hash(minimizer_hash3_rev(&(seq[next_end]),k));   //HASH
            //printf("%d %s %llu, %s %llu\n",next_end, kmers[next_end], invertible_minimizer_hash(minimizer_hash3(kmers[next_end],k)),
            //                             rev_comp, invertible_minimizer_hash(minimizer_hash3(rev_comp,k)));
        }
    }

    *min_l_real = min_position;
//    for(int i = 0; i < kmers_l; i++){
//        printf("%2d h=%22llu, hr = %22llu\n", i, hash_buffer[i], r_hash_buffer[i]);
//    }
//    for(int i = 0; i < *min_l_real; i++){
//        printf("(%lld, %d, %s)\n", (*minimizers)[i].hash, (*minimizers)[i].index, (*minimizers)[i].rev ? "true" : "false");
//    }
//
//    printf("\n\n");

    //FREE BLOK
    delete[] hash_buffer;
    delete[] r_hash_buffer;
    //

    return;
}


void process_sequence(const char* sequence,
                      uint32_t sequence_l,
                      uint32_t w,
                      uint32_t k,
                      std::unordered_multimap<uint64_t, int> &hash_to_index_map,
                      minimizer** ordered_minimizers_addr,
                      uint32_t* ordered_minimizers_l){

    find_minimizers3(sequence, sequence_l, w, k, ordered_minimizers_addr, sequence_l - w - k + 2, ordered_minimizers_l );

    for(uint32_t i = 0; i < *ordered_minimizers_l; i++){
        minimizer tmp = (*ordered_minimizers_addr)[i];
        hash_to_index_map.emplace(tmp.hash, tmp.index);
    }

}



void find_minimizers5
        (const char *seq,
         uint32_t seq_l,
         uint32_t seq_id,
         uint32_t w,
         uint32_t k,
         std::vector<uint64_t >* minimizers,
         std::vector<hashMinPair>* minimizer_hits
        )
{

    const uint32_t kmers_l = seq_l - k + 1;
    //const uint16_t min_l = seq_l - w - k + 2;

    uint64_t *hash_buffer = new uint64_t[w];
    uint64_t *r_hash_buffer = new uint64_t[w];

    for (int i = 0; i < w; i++) {
        hash_buffer[i] = invertible_minimizer_hash(minimizer_hash3(&(seq[i]), k));
        r_hash_buffer[i] = invertible_minimizer_hash(minimizer_hash3_rev(&(seq[i]),k));   //HASH
    }

    uint32_t min_l_pred = seq_l - w - k + 2;
    uint64_t last_min_hash = UINT64_MAX;
    int64_t last_min_position = -1;

    for (uint32_t i = 0; i < min_l_pred; i++) {
        uint64_t u;
        uint64_t v;

        if(last_min_position != -1 && last_min_position >= i && last_min_position < i + w){
            u = hash_buffer[(i + w - 1) % w];
            v = r_hash_buffer[(i + w - 1) % w];

            if(u == v){
                continue;
            }

            else if(u < v && u <= last_min_hash){
                //minimizers->emplace_back((minimizer) {u, (uint32_t)(i + w - 1), false});      //REALLOC
                minimizers->emplace_back(u);
                minimizer_hits->emplace_back((hashMinPair) {u, seq_id, (uint32_t) (i + w - 1), false });
                last_min_position = i + w - 1;
                last_min_hash = u;
            }

            else if(u > v && v <= last_min_hash) {
                //minimizers->emplace_back((minimizer) {v, (uint32_t)(i + w - 1), true});       //REALLOC
                minimizers->emplace_back(v);
                minimizer_hits->emplace_back((hashMinPair) {v, seq_id, (uint32_t) (i + w - 1), true });
                last_min_position = i + w - 1;
                last_min_hash = v;
            }
        }
        else {
            uint64_t m = UINT64_MAX;

            uint32_t *min_positions = new uint32_t[w];
            bool *min_rev = new bool[w];
            uint16_t min_pos_size = 0;
            for (int j = 0; j < w; j++) {
                u = hash_buffer[(i + j) % w];
                v = r_hash_buffer[(i + j) % w];

                if(u == v){
                    continue;
                }

                if(u < m || v <  m){
                    min_pos_size = 0;

                    if(u < v){
                        min_positions[min_pos_size] = i + j;
                        min_rev[min_pos_size] = false;
                        m = u;

                    } else {
                        min_positions[min_pos_size] = i + j;
                        min_rev[min_pos_size] = true;
                        m = v;
                    }

                    min_pos_size++;
                }

                else if(u == m){
                    min_positions[min_pos_size] = i + j;
                    min_rev[min_pos_size++] = false;
                }

                else if(v == m){
                    min_positions[min_pos_size] = i + j;
                    min_rev[min_pos_size++] = true;
                }
            }

            last_min_hash = m;
            last_min_position = min_positions[min_pos_size - 1];

            for(uint32_t j = 0; j < min_pos_size; j++){
                //minimizers->emplace_back((minimizer) {m, (uint32_t)(min_positions[j]), min_rev[j]}); // REALLOC
                minimizers->emplace_back(m);
                minimizer_hits->emplace_back((hashMinPair) {m, seq_id, (uint32_t) (min_positions[j]), min_rev[j] });
            }


            delete[] min_positions;
            delete[] min_rev;
        }

        int next_end = i + w;
        if (next_end < kmers_l) {
            hash_buffer[next_end % w] = invertible_minimizer_hash(minimizer_hash3(&(seq[next_end]),k));//HASH
            r_hash_buffer[next_end % w] = invertible_minimizer_hash(minimizer_hash3_rev(&(seq[next_end]),k));   //HASH
        }
    }

    //FREE BLOK
    delete[] hash_buffer;
    delete[] r_hash_buffer;
    //

    return;
}

void process_sequence2(const char* sequence,
                       uint32_t sequence_l,
                       uint32_t sequence_id,
                       uint32_t w,
                       uint32_t k,
                       std::vector<std::vector<uint64_t>>* ordered_minimizers_addr,
                       std::vector<hashMinPair>* minimizer_hits){

    std::vector<uint64_t> mins;
    find_minimizers5(sequence, sequence_l, sequence_id, w, k, &mins, minimizer_hits);
    mins.shrink_to_fit();
    ordered_minimizers_addr->emplace_back(mins);
}


bool hashMinPair_comparator(hashMinPair a, hashMinPair b){
    if(a.hash < b.hash) return true;
    if(b.hash < a.hash) return false;

    if(a.seq_id < b.seq_id) return true;
    if(b.seq_id < a.seq_id) return  false;

    if(a.index < b.index) return true;
    if(b.index < a.index) return true;
}

bool hashMinPair2_comparator(hashMinPair2 a, hashMinPair2 b){
    return a.index<b.index;
}

void fill_lookup_table(std::vector<hashMinPair>* v, unordered_map<uint64_t, uint64_t>* lookup_table){
    uint64_t last = 0;
    for(uint32_t i = 0, len = v->size(); i < len; i++){
        uint64_t current = (*v)[i].hash;
        if(last == current){
            continue;
        }

        lookup_table->emplace(current, i);
        last = current;
    }
}


void find_minimizers6
        (const char *seq,
         uint32_t seq_l,
         uint32_t seq_id,
         uint32_t w,
         uint32_t k,
         std::vector<uint64_t>& minimizers,
         std::unordered_map<uint64_t, std::vector<hashMinPair2>>& minimizer_hits
        )
{

    const uint32_t kmers_l = seq_l - k + 1;
    //const uint16_t min_l = seq_l - w - k + 2;

    uint64_t *hash_buffer = new uint64_t[w];
    uint64_t *r_hash_buffer = new uint64_t[w];

    for (int i = 0; i < w; i++) {
        hash_buffer[i] = invertible_minimizer_hash(minimizer_hash3(&(seq[i]), k));
        r_hash_buffer[i] = invertible_minimizer_hash(minimizer_hash3_rev(&(seq[i]),k));   //HASH
    }

    uint32_t min_l_pred = seq_l - w - k + 2;
    uint64_t last_min_hash = UINT64_MAX;
    int64_t last_min_position = -1;

    for (uint32_t i = 0; i < min_l_pred; i++) {
        uint64_t u;
        uint64_t v;

        if(last_min_position != -1 && last_min_position >= i && last_min_position < i + w){
            u = hash_buffer[(i + w - 1) % w];
            v = r_hash_buffer[(i + w - 1) % w];

            if(u == v){
                continue;
            }

            else if(u < v && u <= last_min_hash){
                minimizers.emplace_back(u);
                auto it =minimizer_hits.find(u);
                if(it == minimizer_hits.end()){
                    std::vector<hashMinPair2> vec;
                    vec.emplace_back((hashMinPair2) {seq_id, (uint32_t) (i + w - 1), false });
                    minimizer_hits.emplace(u, vec);
                }
                else{
                    it->second.emplace_back((hashMinPair2) {seq_id, (uint32_t) (i + w - 1), false });
                }
                //printf("%llu -> %u, %u, %s\n", u, seq_id, i + w - 1, "True");
                last_min_position = i + w - 1;
                last_min_hash = u;
            }

            else if(u > v && v <= last_min_hash) {
                minimizers.emplace_back(v);
                auto it =minimizer_hits.find(v);
                if(it == minimizer_hits.end()){
                    std::vector<hashMinPair2> vec;
                    vec.emplace_back((hashMinPair2) {seq_id, (uint32_t) (i + w - 1), false });
                    minimizer_hits.emplace(v, vec);
                }
                else{
                    it->second.emplace_back((hashMinPair2) {seq_id, (uint32_t) (i + w - 1), false });
                }
                //printf("%llu -> %u, %u, %s\n", v, seq_id, i + w - 1, "False");
                last_min_position = i + w - 1;
                last_min_hash = v;
            }
        }
        else {
            uint64_t m = UINT64_MAX;

            uint32_t *min_positions = new uint32_t[w];
            bool *min_rev = new bool[w];
            uint16_t min_pos_size = 0;
            for (int j = 0; j < w; j++) {
                u = hash_buffer[(i + j) % w];
                v = r_hash_buffer[(i + j) % w];

                if(u == v){
                    continue;
                }

                if(u < m || v <  m){
                    min_pos_size = 0;

                    if(u < v){
                        min_positions[min_pos_size] = i + j;
                        min_rev[min_pos_size] = false;
                        m = u;

                    } else {
                        min_positions[min_pos_size] = i + j;
                        min_rev[min_pos_size] = true;
                        m = v;
                    }

                    min_pos_size++;
                }

                else if(u == m){
                    min_positions[min_pos_size] = i + j;
                    min_rev[min_pos_size++] = false;
                }

                else if(v == m){
                    min_positions[min_pos_size] = i + j;
                    min_rev[min_pos_size++] = true;
                }
            }

            last_min_hash = m;
            last_min_position = min_positions[min_pos_size - 1];

            for(uint32_t j = 0; j < min_pos_size; j++) {
                minimizers.emplace_back(m);
                auto it = minimizer_hits.find(m);
                if (it == minimizer_hits.end()) {
                    std::vector<hashMinPair2> vec;
                    vec.emplace_back((hashMinPair2) {seq_id, (uint32_t) min_positions[j], min_rev[j]});
                    minimizer_hits.emplace(m, vec);
                } else {
                    it->second.emplace_back((hashMinPair2) {seq_id, (uint32_t) min_positions[j], min_rev[j]});
                }

                //printf("%llu -> %u, %u, %s\n", m, seq_id, min_positions[j], min_rev[j] ? "True" : "False");
            }


            delete[] min_positions;
            delete[] min_rev;
        }

        int next_end = i + w;
        if (next_end < kmers_l) {
            hash_buffer[next_end % w] = invertible_minimizer_hash(minimizer_hash3(&(seq[next_end]),k));//HASH
            r_hash_buffer[next_end % w] = invertible_minimizer_hash(minimizer_hash3_rev(&(seq[next_end]),k));   //HASH
        }
    }

    //FREE BLOK
    delete[] hash_buffer;
    delete[] r_hash_buffer;
    //

    return;
}


void process_sequence3(const char* sequence,
                       uint32_t sequence_l,
                       uint32_t sequence_id,
                       uint32_t w,
                       uint32_t k,
                       std::vector<std::vector<uint64_t>>& ordered_minimizers_addr,
                       std::unordered_map<uint64_t, std::vector<hashMinPair2>>& minimizer_hits){

    std::vector<uint64_t> mins;
    find_minimizers6(sequence, sequence_l, sequence_id, w, k, mins, minimizer_hits);
    mins.shrink_to_fit();
    ordered_minimizers_addr.emplace_back(mins);
}

void sort_by_indices(std::unordered_map<uint64_t, std::vector<hashMinPair2>>& minimizer_hits) {
    auto it = minimizer_hits.begin();

    std::shared_ptr<thread_pool::ThreadPool> thread_pool =
            thread_pool::createThreadPool();
    std::vector<std::future<void>> thread_futures;

    //sortiranje vektora u mapi i concurrency??
    while (it != minimizer_hits.end()) {
        auto a = it-> second.begin();
        thread_futures.emplace_back(thread_pool->submit_task(
                std::sort, it->second.begin(), it->second.end(), hashMinPair2_comparator));
        it++;
    }

    for (auto &it: thread_futures) {
        it.wait();
    }
}


void find_minimizers7
        (const char *seq,
         uint32_t seq_l,
         uint32_t seq_id,
         uint32_t w,
         uint32_t k,
         std::vector<minimizer>& minimizers
        )
{

    const uint32_t kmers_l = seq_l - k + 1;
    //const uint16_t min_l = seq_l - w - k + 2;

    uint64_t *hash_buffer = new uint64_t[w];
    uint64_t *r_hash_buffer = new uint64_t[w];

    for (int i = 0; i < w; i++) {
        hash_buffer[i] = invertible_minimizer_hash(minimizer_hash3(&(seq[i]), k));
        r_hash_buffer[i] = invertible_minimizer_hash(minimizer_hash3_rev(&(seq[i]),k));   //HASH
    }

    uint32_t min_l_pred = seq_l - w - k + 2;
    uint64_t last_min_hash = UINT64_MAX;
    int64_t last_min_position = -1;

    for (uint32_t i = 0; i < min_l_pred; i++) {
        uint64_t u;
        uint64_t v;

        if(last_min_position != -1 && last_min_position >= i && last_min_position < i + w){
            u = hash_buffer[(i + w - 1) % w];
            v = r_hash_buffer[(i + w - 1) % w];

            if(u == v){
                continue;
            }

            else if(u < v && u <= last_min_hash){
                minimizers.emplace_back((minimizer) {u, (uint32_t) (i + w - 1), false });
                //printf("%llu -> %u, %u, %s\n", u, seq_id, i + w - 1, "True");
                last_min_position = i + w - 1;
                last_min_hash = u;
            }

            else if(u > v && v <= last_min_hash) {
                minimizers.emplace_back((minimizer) {v, (uint32_t) (i + w - 1), false });
                //printf("%llu -> %u, %u, %s\n", v, seq_id, i + w - 1, "False");
                last_min_position = i + w - 1;
                last_min_hash = v;
            }
        }
        else {
            uint64_t m = UINT64_MAX;

            uint32_t *min_positions = new uint32_t[w];
            bool *min_rev = new bool[w];
            uint16_t min_pos_size = 0;
            for (int j = 0; j < w; j++) {
                u = hash_buffer[(i + j) % w];
                v = r_hash_buffer[(i + j) % w];

                if(u == v){
                    continue;
                }

                if(u < m || v <  m){
                    min_pos_size = 0;

                    if(u < v){
                        min_positions[min_pos_size] = i + j;
                        min_rev[min_pos_size] = false;
                        m = u;

                    } else {
                        min_positions[min_pos_size] = i + j;
                        min_rev[min_pos_size] = true;
                        m = v;
                    }

                    min_pos_size++;
                }

                else if(u == m){
                    min_positions[min_pos_size] = i + j;
                    min_rev[min_pos_size++] = false;
                }

                else if(v == m){
                    min_positions[min_pos_size] = i + j;
                    min_rev[min_pos_size++] = true;
                }
            }

            last_min_hash = m;
            last_min_position = min_positions[min_pos_size - 1];

            for(uint32_t j = 0; j < min_pos_size; j++) {
                minimizers.emplace_back((minimizer) {m, (uint32_t) min_positions[j], min_rev[j]});

                //printf("%llu -> %u, %u, %s\n", m, seq_id, min_positions[j], min_rev[j] ? "True" : "False");
            }


            delete[] min_positions;
            delete[] min_rev;
        }

        int next_end = i + w;
        if (next_end < kmers_l) {
            hash_buffer[next_end % w] = invertible_minimizer_hash(minimizer_hash3(&(seq[next_end]),k));//HASH
            r_hash_buffer[next_end % w] = invertible_minimizer_hash(minimizer_hash3_rev(&(seq[next_end]),k));   //HASH
        }
    }

    //FREE BLOK
    delete[] hash_buffer;
    delete[] r_hash_buffer;
    //

    return;
}


void process_sequence4(const char* sequence,
                       uint32_t sequence_l,
                       uint32_t sequence_id,
                       uint32_t w,
                       uint32_t k,
                       std::vector<std::vector<minimizer>>& ordered_minimizers_addr){

    std::vector<minimizer> mins;
    find_minimizers7(sequence, sequence_l, sequence_id, w, k, mins);
    mins.shrink_to_fit();
    ordered_minimizers_addr[sequence_id] = mins;
}

void fill_lookup_table(std::vector<std::vector<minimizer>> minimizers, std::unordered_map<uint64_t, vector<hashMinPair2>>& map){
    for(uint32_t i = 0, len = minimizers.size(); i < len; i++){
        vector<minimizer> tmp = minimizers[i];
        for(uint32_t j = 0, len2 = tmp.size(); j < len2; j++){
            minimizer tmp2 = tmp[j];
            auto it = map.find(tmp2.hash);
            if (it == map.end()) {
                std::vector<hashMinPair2> vec;
                vec.emplace_back((hashMinPair2) {i, tmp2.index, tmp2.rev});
                map.emplace(tmp2.hash, vec);
            } else {
                it->second.emplace_back((hashMinPair2) {i, tmp2.index, tmp2.rev});
            }
        }
    }
}

void fill_lookup_table_nogo_minimizers(std::vector<std::vector<minimizer>>& minimizers,
                                       std::unordered_map<uint64_t,
                                               vector<hashMinPair2>>& map,
                                       std::vector<uint64_t>& no_gos,
                                       double threshold){

    std::unordered_map<uint64_t ,uint32_t> min_occurences(map.size());
    uint64_t num_of_minimizers = 0;

    for(uint32_t i = 0, len = minimizers.size(); i < len; i++){
        vector<minimizer> tmp = minimizers[i];
        num_of_minimizers += tmp.size();

        for(uint32_t j = 0, len2 = tmp.size(); j < len2; j++){
            minimizer tmp2 = tmp[j];
            auto it = map.find(tmp2.hash);
            if (it == map.end()) {
                std::vector<hashMinPair2> vec;
                vec.emplace_back((hashMinPair2) {i, tmp2.index, tmp2.rev});
                map.emplace(tmp2.hash, vec);
            } else {
                it->second.emplace_back((hashMinPair2) {i, tmp2.index, tmp2.rev});
            }

            min_occurences[tmp2.hash]++;
        }
    }

    std::vector<std::pair<uint64_t, uint32_t>> min_occur(min_occurences.begin(), min_occurences.end());
    sort(min_occur.begin(), min_occur.end(), occurences_comparator);
    double acc = 0;
    for(int i = 0, len = min_occur.size(); i < len; i++){
        acc += min_occur[i].second / (double) num_of_minimizers;

        if(acc > threshold){
            break;
        }
        no_gos.emplace_back(min_occur[i].first);
    }


}

bool occurences_comparator(std::pair<uint64_t,uint32_t>& a, std::pair<uint64_t,uint32_t>& b){
    return a.second > b.second;
}