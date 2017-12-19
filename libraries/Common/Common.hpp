
#ifndef Common
#define Common

#include <vector>
#include <string>
#include <functional>
#include <unordered_map>
#include <CustomTypes.h>

using namespace std;

typedef tuple<string, int, string> triplet;

typedef tuple<uint64_t, int, int> triple;

std::vector<std::string> find_kmer(int k, std::string s);

uint64_t minimizer_hash(string s);

vector<triplet> find_minimizers(int w, int k, string s);

vector<tuple<uint64_t, int, int>> find_minimizers2(int w, int k, string s);

size_t find_hash_value(char c);

std::unordered_multimap<uint64_t, tuple<string, int, int>, function<size_t( uint64_t)>> indexSequences(vector<string> sequences, int w, int k);

std::unordered_multimap<uint64_t, hashEntry, function<size_t( uint64_t)>> indexSequence(string sequence, int w, int k);

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

const int hash_width = 20;

const int offset_width = 6;
typedef tuple<uint64_t, tuple<string, int, int>> hashMinPair3;

std::vector<hashMinPair3> indexTable(vector<string> sequences, int w, int k);

typedef tuple<string, int, int, int> minimizer_hit;

bool hit_comparator(const minimizer_hit a,const minimizer_hit b);

vector<mapInfo> map_minimizers(unordered_multimap<uint64_t, tuple<string, int, int>, function<size_t(uint64_t)>> lookup_table, string query_sequence, int w, int k, int epsilon);

//new optimized functions


char complement(char c);
void find_reverse_complement(const char *seq, uint32_t seq_l, char** rev_comp);
void destroy_reverse_complement(char** rev_comp);
void find_kmers(const char *seq, uint32_t k, char ***kmers, uint32_t kmers_l);
void destroy_kmers(char*** kmers, uint32_t kmers_l);
uint64_t minimizer_hash3(const char* seq, uint32_t seq_l);
uint64_t minimizer_hash3_rev(const char* seq, uint32_t seq_l);
uint64_t invertible_minimizer_hash(uint64_t x);
uint64_t invertible_minimizer_hash_inverse(uint64_t x);
void find_minimizers3(const char *seq, uint32_t seq_l, uint32_t w, uint32_t k, minimizer** minimizers, uint32_t min_l_pred, uint32_t* min_l_real);
void find_minimizers4(const char *seq, uint32_t seq_l, uint32_t w, uint32_t k, minimizer** minimizers, uint32_t min_l_pred, uint32_t* min_l_real);
void find_minimizers5
        (const char *seq,
         uint32_t seq_l,
         uint32_t seq_id,
         uint32_t w,
         uint32_t k,
         std::vector<uint64_t >* minimizers,
         std::vector<hashMinPair>* minimizer_hits
        );
void process_sequence2(const char* sequence,
                       uint32_t sequence_l,
                       uint32_t sequence_id,
                       uint32_t w,
                       uint32_t k,
                       std::vector<std::vector<uint64_t >>* ordered_minimizers_addr,
                       std::vector<hashMinPair>* minimizer_hits);
bool hashMinPair_comparator(hashMinPair a, hashMinPair b);
void fill_lookup_table(std::vector<hashMinPair>* v, unordered_map<uint64_t, uint64_t>* lookup_table);

void process_sequence3(const char* sequence,
                       uint32_t sequence_l,
                       uint32_t sequence_id,
                       uint32_t w,
                       uint32_t k,
                       std::vector<std::vector<uint64_t>>& ordered_minimizers_addr,
                       std::unordered_map<uint64_t, std::vector<hashMinPair2>>& minimizer_hits);
void find_minimizers6
        (const char *seq,
         uint32_t seq_l,
         uint32_t seq_id,
         uint32_t w,
         uint32_t k,
         std::vector<uint64_t >& minimizers,
         std::unordered_map<uint64_t, std::vector<hashMinPair2>>& minimizer_hits
        );

void sort_by_indices(std::unordered_map<uint64_t, std::vector<hashMinPair2>>& minimizer_hits);

#endif
