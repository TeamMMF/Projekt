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

std::vector<std::pair<int, bool>> find_overlaps_by_LIS_parallel(int  query_id,
                                                                std::vector<minimizer>& minimizers,
                                                                std::unordered_map<uint64_t, std::vector<hashMinPair2>>&  minimizers_for_hash,
                                                                int lis_threshold,
                                                                std::vector<uint64_t>& nogos);

std::vector<std::pair<int, bool>> find_overlaps_by_LIS_parallel(int  query_id,
                                                           std::vector<minimizer>& minimizers,
                                                           std::unordered_map<uint64_t, std::vector<hashMinPair2>>&  minimizers_for_hash,
                                                           int lis_threshold,
                                                           std::unordered_map<uint64_t, uint32_t >& occurrences);

#endif