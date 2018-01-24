
#ifndef Common
#define Common

#include <vector>
#include <string>
#include <functional>
#include <unordered_map>
#include <CustomTypes.h>

using namespace std;


size_t find_hash_value(char c);

/*
 * sequence - sekvenca za koju trazimo minimizere
 * w, k -klasika
 * hash_to_index_map_addr - adresa na koju se treba staviti mapa hashMinimizera -> index
 * ordered_minimizers_addr - adresa na koju se treba staviti sortirano polje parova (minimizer_index, minimizer_hash)
 */
void process_sequence(const char* sequence,
                      uint32_t sequence_l,
                      uint32_t w,
                      uint32_t k,
                      std::unordered_multimap<uint64_t, int> &hash_to_index_map,
                      minimizer** ordered_minimizers_addr,
                      uint32_t* ordered_minimizers_l);

typedef tuple<string, int, int, int> minimizer_hit;

bool hit_comparator(const minimizer_hit a,const minimizer_hit b);




char complement(char c);

uint64_t minimizer_hash3(const char* seq, uint32_t seq_l);
uint64_t minimizer_hash3_rev(const char* seq, uint32_t seq_l);
uint64_t invertible_minimizer_hash(uint64_t x);
uint64_t invertible_minimizer_hash_inverse(uint64_t x);

void process_sequence3(const char* sequence,
                       uint32_t sequence_l,
                       uint32_t sequence_id,
                       uint32_t w,
                       uint32_t k,
                       std::vector<std::vector<uint64_t>>& ordered_minimizers_addr,
                       std::unordered_map<uint64_t, std::vector<hashMinPair2>>& minimizer_hits);


void sort_by_indices(std::unordered_map<uint64_t, std::vector<hashMinPair3>>& minimizer_hits);

void process_sequence4(const char* sequence,
                       uint32_t sequence_l,
                       uint32_t sequence_id,
                       int32_t w,
                       uint32_t k,
                       std::vector<std::vector<minim>>& ordered_minimizers_addr);

uint32_t process_sequence4_id(const char* sequence,
                              uint32_t sequence_l,
                              uint32_t sequence_id,
                              int32_t w,
                              uint32_t k,
                              std::vector<std::vector<minim>>& ordered_minimizers_addr);
void find_minimizers7
        (const char *seq,
         uint32_t seq_l,
         uint32_t seq_id,
         int32_t w,
         uint32_t k,
         std::vector<minim>& minimizers);


void process_sequence_all(const char* sequence,
                          uint32_t sequence_l,
                          uint32_t sequence_id,
                          int32_t w,
                          uint32_t k,
                          std::vector<std::vector<minim>>& ordered_minimizers_addr,
                          std::unordered_map<uint64_t, std::vector<hashMinPair3>>& minimizer_hits,
                          std::unordered_map<uint64_t,uint32_t>& occurrences);


bool occurences_comparator(const std::pair<uint64_t,uint32_t>& a, const std::pair<uint64_t,uint32_t>& b);

void fill_lookup_table_nogo_minimizers(std::vector<std::vector<minim>>& minimizers, std::unordered_map<uint64_t, vector<hashMinPair3>>& map,
                                       std::vector<uint64_t>& no_gos, double threshold);
void sort_wrap(vector<hashMinPair2>::iterator  begin, vector<hashMinPair2>::iterator end);

bool hashMinPair3_comparator(hashMinPair3 a, hashMinPair3 b);

void find_maximizers7
        (const char *seq,
         uint32_t seq_l,
         uint32_t seq_id,
         int32_t w,
         uint32_t k,
         std::vector<minim>& minimizers
        );

void find_maximizers_full
        (const char *seq,
         uint32_t seq_l,
         uint32_t seq_id,
         int32_t w,
         uint32_t k,
         std::vector<minim>& minimizers,
         std::unordered_map<uint64_t, std::vector<hashMinPair3>>& minimizer_hits,
         std::unordered_map<uint64_t ,uint32_t >& occurrences
        );

void process_sequence4_max(const char* sequence,
                           uint32_t sequence_l,
                           uint32_t sequence_id,
                           int32_t w,
                           uint32_t k,
                           std::vector<std::vector<minim>>& ordered_minimizers_addr);

uint64_t minimizer_hash4(const char* seq, int32_t index, uint64_t* last_hash, uint32_t power, uint64_t* first_nucleotide_value);
uint64_t minimizer_hash4_rev(const char* seq, int32_t index, uint64_t* last_hash, uint32_t power, uint64_t* first_nucleotide_value);
#endif
