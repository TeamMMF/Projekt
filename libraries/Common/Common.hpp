
#ifndef Common
#define Common

#include <vector>
#include <string>
#include <functional>
#include <unordered_map>
#include <CustomTypes.h>

using namespace std;


size_t find_hash_value(char c);
typedef tuple<string, int, int, int> minimizer_hit;

bool hit_comparator(const minimizer_hit a,const minimizer_hit b);




char complement(char c);

uint64_t minimizer_hash3(const char* seq, uint32_t seq_l);
uint64_t minimizer_hash3_rev(const char* seq, uint32_t seq_l);
uint64_t invertible_minimizer_hash(uint64_t x);
uint64_t invertible_minimizer_hash_inverse(uint64_t x);

void sort_by_indices(std::unordered_map<uint64_t, std::vector<hashMinPair3>>& minimizer_hits);

uint32_t process_sequence4_id(const char* sequence,
                              uint32_t sequence_l,
                              uint32_t sequence_id,
                              int32_t w,
                              uint32_t k,
                              std::vector<std::vector<minim>>& ordered_minimizers_addr);
void find_minimizers7
        (const char *seq,
         uint32_t seq_l,
         int32_t w,
         uint32_t k,
         std::vector<minim>& minimizers);

bool hashMinPair3_comparator(hashMinPair3 a, hashMinPair3 b);

uint64_t minimizer_hash4(const char* seq, int32_t index, uint64_t* last_hash, uint32_t power, uint64_t* first_nucleotide_value);
uint64_t minimizer_hash4_rev(const char* seq, int32_t index, uint64_t* last_hash, uint32_t power, uint64_t* first_nucleotide_value);
#endif
