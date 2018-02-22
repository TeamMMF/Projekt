//
// Created by matep on 03/12/2017.
//
#ifndef Dynamic
#define Dynamic
#define IGNORE_THRESHOLD 34

#include <tuple>
#include <string>
#include <vector>
#include <cstdint>
#include "CustomTypes.h"
#include <unordered_map>


std::vector<std::pair<int, bool>> find_overlaps_by_LIS(int query_id,
                                                       std::vector<minim> &minimizers,
                                                       std::unordered_map<uint64_t, std::vector<hashMinPair3>> &lookup_map,
                                                       int lis_threshold,
                                                       std::unordered_map<uint64_t, uint32_t> &occurrences);

#endif