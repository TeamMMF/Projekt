//
// Created by matep on 03/12/2017.
//
#ifndef Dynamic
#define Dynamic
#include <tuple>
#include <string>
#include <vector>
#include <cstdint>
#include "CustomTypes.h"
#include <unordered_map>

void LCS(const char *s1, const char *s2, int s1_l, int s2_l);

void print_matrix(const char *s1, const char *s2, uint8_t *matrix, int r, int c);

bool check_k_substring_match(const char* s1, const char* s2, int i, int j, int k);

void LCS_k(const char* s1, const char* s2, int s1_l, int s2_l, int k);

std::string traceback_LCS(uint8_t *matrix, int r, int c, const char* s1);

typedef std::tuple<int, int, bool> matchPair;

std::vector<matchPair> generate_match_pairs(std::string s1, std::string s2, int k);

bool matchPair_comparator(const matchPair a, const matchPair b);

uint64_t max_between_indexes(std::vector<uint64_t> array, int start, int end);

bool check_precedence(matchPair p, matchPair g);

uint64_t LCS_kpp(std::string s1, std::string s2, int k);

size_t tuple_hash(std::tuple<int, int, bool> x);

int lis(const std::vector<int> &v);

int compare_with_lis(minimizer* seq1_mins_sorted,
                                     int seq1_mins_size,
                                     std::unordered_multimap<uint64_t, int> &seq2_hash_to_index,
                                     minimizer* seq2_mins_sorted,bool same_strand);

std::pair<int,char> compare_with_lis(minimizer* seq1_mins_sorted,
                     int seq1_mins_size,
                     std::unordered_multimap<uint64_t, int> &seq2_hash_to_index,
                     minimizer* seq2_mins_sorted);

std::vector<std::pair<int, bool>> find_overlaps_by_LIS(int  query_id,
                                                  std::vector<uint64_t >& minimizer_hashes,
                                                  std::unordered_map<uint64_t, std::vector<hashMinPair2>>& minimizers_for_hash,
                                                  int lis_threshold);

std::vector<std::pair<int, bool>> find_overlaps_by_LIS_parallel(int  query_id,
                                             std::vector<minimizer>& minimizers,
                                             std::unordered_map<uint64_t, std::vector<hashMinPair2>>&  minimizers_for_hash,
                                             int lis_threshold);

#endif