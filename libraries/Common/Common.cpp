

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
/**
 * Returns the associated hash value of a given nucleotide
 * @param c nucleotide representation
 * @return nucleotide hash value
 */
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
            return 0;
    }
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
 * Naive implementation of the minimizer hashing function. The kmer hash is calculated based on the position of each
 * nucleotide in the sequence as well as its hash value.
 * @param seq  pointer to the start of the nucleotide char sequence
 * @param seq_l length of the kmer
 * @return 64-bit hash value
 */
uint64_t minimizer_hash3(const char* seq, uint32_t seq_l)
{
    uint64_t hash = 0;

    for (int i = 0; i < seq_l; i++) {
        hash += find_hash_value(seq[i]) * pow(4, seq_l - i - 1);
    }

    return hash;
}

/**
 * Naive implementation of the minimizer hashing function. The kmer hash is calculated based on the position of each
 * nucleotide in the reversed sequence as well as its hash value.
 * @param seq  pointer to the start of the nucleotide char sequence
 * @param seq_l length of the kmer
 * @return 64-bit hash value
 */
uint64_t minimizer_hash3_rev(const char* seq, uint32_t seq_l) {
    uint64_t hash = 0;

    for (int i = seq_l - 1; i >= 0; i--) {
        hash += find_hash_value(complement(seq[i])) * pow(4, i);
    }

    return hash;
}


/**
 * Compares two hashMinPair3 structs based on the index.
 * @param a first hashMinPair3 struct
 * @param b second hashMinPair3 struct
 * @return true if absolute value of the first struct index is less than the absolute value of the second struct index
 */
bool hashMinPair3_comparator(hashMinPair3 a, hashMinPair3 b){
    return abs(a.index)<abs(b.index);
}


/**
 * Sorts each vector using the hahsMinPair3_comparator of the given unoredered_map
 * @param minimizer_hits reference to unordered_map
 */
void sort_by_indices(std::unordered_map<uint64_t, std::vector<hashMinPair3>>& minimizer_hits) {
    auto it = minimizer_hits.begin();
    while (it != minimizer_hits.end()) {
        sort(it->second.begin(), it->second.end(),hashMinPair3_comparator);
        it++;
    }
}

/**
 * Find all minimizers for a given strand and stores them as a minim struct in a vector.
 * @param seq sequence
 * @param seq_l sequence length
 * @param w window size
 * @param k kmer length
 * @param minimizers  reference to vector where all minimizers will be stored
 */
void find_minimizers7
        (const char *seq,
         uint32_t seq_l,
         int32_t w,
         uint32_t k,
         std::vector<minim>& minimizers
        )
{

    if(seq_l < w + k){
        return;
    }

    const uint32_t kmers_l = seq_l - k + 1;

    uint64_t *hash_buffer = new uint64_t[w];
    uint64_t *r_hash_buffer = new uint64_t[w];
    uint64_t first_nuc_val = find_hash_value(seq[0]);
    uint64_t first_nuc_val_r = find_hash_value(complement(seq[0]));
    uint32_t power = k - 1;
    uint64_t prev_hash = minimizer_hash3(seq, k);
    uint64_t prev_hash_r = minimizer_hash3_rev(seq,k);

    hash_buffer[0] = invertible_minimizer_hash(prev_hash);
    r_hash_buffer[0] = invertible_minimizer_hash(prev_hash_r);

    for (uint32_t i = 1; i < w; i++) {
        hash_buffer[i] = invertible_minimizer_hash(minimizer_hash4(seq, i, &prev_hash, power, &first_nuc_val));
        r_hash_buffer[i] = invertible_minimizer_hash(minimizer_hash4_rev(seq, i, &prev_hash_r, power, &first_nuc_val_r));   //HASH
    }

    uint32_t min_l_pred = seq_l - w - k + 2;
    uint64_t last_min_hash = UINT64_MAX;
    int64_t last_min_position = -1;

    for (int32_t i = 0; i < min_l_pred; i++) {
        uint64_t u;
        uint64_t v;

        if(last_min_position != -1 && abs(last_min_position) >= i){
            u = hash_buffer[(i + w - 1) % w];
            v = r_hash_buffer[(i + w - 1) % w];

            if(u == v){
                continue;
            }

            else if(u < v && u <= last_min_hash){
                minimizers.emplace_back((minim) {u, (i + w - 1)});
                last_min_position = i + w - 1;
                last_min_hash = u;
            }

            else if(u > v && v <= last_min_hash) {
                minimizers.emplace_back((minim) {v, -(i + w - 1)});
                last_min_position = -(i + w - 1);
                last_min_hash = v;
            }
        }
        else {
            uint64_t m = UINT64_MAX;

            int32_t *min_positions = new int32_t[w];
            uint16_t min_pos_size = 0;
            for (int j = 0; j < w; j++) {
                u = hash_buffer[(i + j) % w];
                v = r_hash_buffer[(i + j) % w];

                if(u == v){
                    continue;
                }

                if(u < m || v <  m){
                    if(u < v){
                        min_positions[0] = i + j;
                        m = u;

                    } else {
                        min_positions[0] = (i + j);
                        m = v;
                    }

                    min_pos_size = 1;
                }

                else if(u == m){
                    min_positions[min_pos_size++] = i + j;
                }

                else if(v == m){
                    min_positions[min_pos_size++] = -(i + j);
                }
            }

            last_min_hash = m;
            last_min_position = min_positions[min_pos_size - 1];

            for(uint32_t j = 0; j < min_pos_size; j++) {
                minimizers.emplace_back((minim) {m,  min_positions[j]}); //TREBA NEGATIVNO
            }


            delete[] min_positions;
        }

        int next_end = i + w;
        if (next_end < kmers_l) {
            hash_buffer[next_end % w] = invertible_minimizer_hash(minimizer_hash4(seq, next_end, &prev_hash, power, &first_nuc_val));//HASH
            r_hash_buffer[next_end % w] = invertible_minimizer_hash(minimizer_hash4_rev(seq, next_end, &prev_hash_r, power, &first_nuc_val_r));   //HASH
        }
    }

    delete[] hash_buffer;
    delete[] r_hash_buffer;

    return;
}

/**
 * Wrapper function for finding all minimizers of a sequence and embedding the resulting vector into a vector of vectors
 * which keeps track of all minimizer vectors for all sequnces.
 * @param sequence nucleotide sequence from which minimizers should be extracted
 * @param sequence_l sequence length
 * @param sequence_id sequence identifier
 * @param w window size
 * @param k kmers size
 * @param ordered_minimizers_addr reference to vector where the reseulting minimizer vector should be stored
 * @return sequence id
 */
uint32_t process_sequence4_id(const char* sequence,
                           uint32_t sequence_l,
                           uint32_t sequence_id,
                           int32_t w,
                           uint32_t k,
                           std::vector<std::vector<minim>>& ordered_minimizers_addr){

    std::vector<minim> mins;
    find_minimizers7(sequence, sequence_l, w, k, mins);
    mins.shrink_to_fit();
    ordered_minimizers_addr[sequence_id] = mins;
    return sequence_id;
}



/**
 * Calculates the minimizer hash of a kmer based on the value of the kmer before it and the last current nucleotide
 * @param seq sequence
 * @param index start index of the kmer in the sequence
 * @param last_hash pointer to last recorded kmer hash valie
 * @param power power to which values should be raised
 * @param first_nucleotide_value poitner to the first value of the nucleotide
 * @return new hash value
 */
uint64_t minimizer_hash4(const char* seq, int32_t index, uint64_t*  last_hash, uint32_t power, uint64_t* first_nucleotide_value)
{
    uint64_t  hash = 0;
    uint64_t remove = *first_nucleotide_value * (1 << (power*2));
    *first_nucleotide_value = find_hash_value(seq[index]);
    hash = (*last_hash - remove) * 4 + find_hash_value(seq[index + power]);
    *last_hash = hash;
    return hash;
}

/**
 * Calculates the reverse minimizer hash of a kmer based on the value of the kmer before it and the last current nucleotide
 * @param seq sequence
 * @param index start index of the kmer in the sequence
 * @param last_hash pointer to last recorded kmer hash valie
 * @param power power to which values should be raised
 * @param first_nucleotide_value poitner to the first value of the nucleotide
 * @return new hash value
 */
uint64_t minimizer_hash4_rev(const char* seq, int32_t index, uint64_t*  last_hash, uint32_t power, uint64_t* first_nucleotide_value)
{
    uint64_t  hash = 0;
    uint64_t remove = *first_nucleotide_value ;
    *first_nucleotide_value = find_hash_value(complement(seq[index]));
    hash = (*last_hash - remove) / 4 + (find_hash_value(complement(seq[index + power])) * (1 << (power*2)));
    *last_hash = hash;
    return hash;
}
