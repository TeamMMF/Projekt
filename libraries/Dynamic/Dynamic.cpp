//
// Created by matep on 03/12/2017.
//
#include "Dynamic.h"
#include "Common.hpp"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <climits>
#include <CustomTypes.h>



using namespace std;



int compare_with_lis(minimizer* seq1_mins_sorted,
                     int seq1_mins_size,
                     unordered_multimap<uint64_t, int> &seq2_hash_to_index,
                     minimizer* seq2_mins_sorted, bool same_strand){
    vector<int> lis_arr;

    for(int i = 0; i < seq1_mins_size; i++){
        auto range = seq2_hash_to_index.equal_range(seq1_mins_sorted[i].hash);
        if(range.first._M_cur == NULL)
            continue;
        int min_diff = INT_MAX;
        int final = -1;
        for (auto it = range.first; it != range.second; ++it) {
            int cur_diff = abs(it->second - seq1_mins_sorted[i].index);
            if((seq2_mins_sorted[it -> second].rev == seq1_mins_sorted[i].rev && same_strand
                || seq2_mins_sorted[it -> second].rev != seq1_mins_sorted[i].rev && !same_strand)
               && cur_diff < min_diff){
                final = it -> second;
                min_diff = cur_diff;
            }
        }
        if(final!=-1) {
            lis_arr.push_back(final);
        }
    }

    if(!same_strand)
        reverse(lis_arr.begin(),lis_arr.end());

    return lis(lis_arr);

}



pair<int,char> compare_with_lis(minimizer* seq1_mins_sorted,
                                int seq1_mins_size,
                                unordered_multimap<uint64_t, int> &seq2_hash_to_index,
                                minimizer* seq2_mins_sorted){

    vector<int> lis_arr_same;
    vector<int> lis_arr_diff;

    //fprintf(stdout, "SEDMI SPRAT\n");
    for(int i = 0; i < seq1_mins_size; i++){

        auto range = seq2_hash_to_index.equal_range(seq1_mins_sorted[i].hash);
        if(range.first._M_cur == nullptr)
            continue;

        int min_diff_same = INT_MAX;
        int min_diff_diff = INT_MAX;
        int final_same =range.first->second; // zasto ovo radi bolje ovako nego kad je tu -1??
        int final_diff =range.first->second;

        for (auto it = range.first; it != range.second; ++it) {
            int cur_diff = (int) abs(it->second - seq1_mins_sorted[i].index);

            if(seq2_mins_sorted[it -> second].rev == seq1_mins_sorted[i].rev && cur_diff < min_diff_same){
                final_same = it -> second;
                min_diff_same = cur_diff;

            }else if (seq2_mins_sorted[it -> second].rev != seq1_mins_sorted[i].rev && cur_diff < min_diff_diff){
                final_diff = it -> second;
                min_diff_diff = cur_diff;
            }
        }

        if(final_same!=-1) {
            lis_arr_same.push_back(final_same);
        }
        if(final_diff!=-1) {
            lis_arr_diff.push_back(final_diff);
        }
    }

    reverse(lis_arr_diff.begin(),lis_arr_diff.end());

    int lis_same = lis (lis_arr_same);
    int lis_diff = lis (lis_arr_diff);

    if(lis_diff > lis_same){
        return make_pair(lis_diff,'-');
    }else{
        return make_pair(lis_same,'+');
    }

}

vector<pair<int, bool>> find_overlaps_by_LIS(int  query_id,
                                             vector<uint64_t >& minimizer_hashes,
                                             unordered_map<uint64_t,
                                                     vector<hashMinPair2>>&minimizers_for_hash,
                                             int lis_threshold){
    unordered_map<uint64_t, vector<int>> same_strand;
    unordered_map<uint64_t, vector<int>> different_strand;

    for(auto h : minimizer_hashes){
        auto matches = minimizers_for_hash.find(h);
        if(matches == minimizers_for_hash.end())
            continue;
        bool curr_rev;
        for(auto match : matches->second){
            if(match.seq_id == query_id){
                curr_rev = match.rev;
            }
            break;
        }
        for(auto match : matches->second){
            if(match.seq_id <= query_id){
                continue;
            }
            if(match.rev ^ curr_rev){
                different_strand[match.seq_id].push_back(match.index);
            }else{
                same_strand[match.seq_id].push_back(match.index);
            }
        }
    }

    vector<pair<int, bool>> overlaps;
    for(auto &entry : same_strand){
        if(lis(entry.second)>=lis_threshold){
            overlaps.emplace_back(entry.first,true);
        }
    }

    for(auto &entry : different_strand){
        reverse(entry.second.begin(),entry.second.end());
        if(lis(entry.second)>=lis_threshold){
            overlaps.emplace_back(make_pair(entry.first,false));
        }
    }

    return overlaps;
};

vector<pair<int, bool>> find_overlaps_by_LIS_parallel(int  query_id,
                                             vector<minimizer>& minimizers,
                                             unordered_map<uint64_t,
                                                     vector<hashMinPair2>>&  minimizers_for_hash,
                                             int lis_threshold){
    unordered_map<uint64_t, vector<int>> same_strand;
    unordered_map<uint64_t, vector<int>> different_strand;
    for(auto min : minimizers){
        auto matches = minimizers_for_hash.find(min.hash);
        if(matches == minimizers_for_hash.end())
            continue;
        for(auto match : matches->second){
            if(match.seq_id <= query_id){
                continue;
            }
            if(match.rev ^ min.rev){
                different_strand[match.seq_id].push_back(-match.index);
            }else{
                same_strand[match.seq_id].push_back(match.index);
            }
        }
    }
    vector<pair<int, bool>> overlaps;
    for(auto &entry : same_strand){
        if(lis(entry.second)>=lis_threshold){
            overlaps.emplace_back(entry.first,true);
        }
    }

    for(auto &entry : different_strand){
        if(lis(entry.second)>=lis_threshold){
            overlaps.emplace_back(make_pair(entry.first,false));
        }
    }
    return overlaps;
};


vector<pair<int, bool>> find_overlaps_by_LIS_parallel(int  query_id,
                                                      vector<minim>& minimizers,
                                                      unordered_map<uint64_t, vector<hashMinPair3>>&  minimizers_for_hash,
                                                      int lis_threshold,
                                                      vector<uint64_t>& nogos){
    unordered_map<uint64_t, vector<int>> same_strand;
    unordered_map<uint64_t, vector<int>> different_strand;
    for(auto min : minimizers){
        if(binary_search(nogos.begin(),nogos.end(),min.hash))
            continue;
        auto matches = minimizers_for_hash.find(min.hash);
        if(matches == minimizers_for_hash.end())
            continue;
        for(auto match : matches->second){
            if(match.seq_id <= query_id){
                continue;
            }
            int final = abs(match.index);
            if(match.index * min.index < 0 ){
                different_strand[match.seq_id].push_back(-final);
            }else{
                same_strand[match.seq_id].push_back(final);
            }
        }
    }
    vector<pair<int, bool>> overlaps;
    for(auto &entry : same_strand){
        if(lis(entry.second)>=lis_threshold){
            overlaps.emplace_back(entry.first,true);
        }
    }

    for(auto &entry : different_strand){
        if(lis(entry.second)>=lis_threshold){
            overlaps.emplace_back(make_pair(entry.first,false));
        }
    }
    return overlaps;
};


vector<pair<int, bool>> find_overlaps_by_LIS_parallel(int  query_id,
                                                      vector<minim>& minimizers,
                                                      unordered_map<uint64_t, vector<hashMinPair3>>&  minimizers_for_hash,
                                                      int lis_threshold,
                                                      unordered_map<uint64_t, uint32_t >& occurrences){
    unordered_map<uint64_t, vector<int>> same_strand;
    unordered_map<uint64_t, vector<int>> different_strand;
    for(auto min : minimizers){
        if(occurrences[min.hash]>34)
            continue;
        auto matches = minimizers_for_hash.find(min.hash);
        if(matches == minimizers_for_hash.end())
            continue;
        for(auto match : matches->second){
            if(match.seq_id > query_id){
                break;
            }
            int final = abs(match.index);
            if(match.index * min.index < 0){
                different_strand[match.seq_id].push_back(-final);
            }else{
                same_strand[match.seq_id].push_back(final);
            }
        }
    }
    vector<pair<int, bool>> overlaps;
    for(auto &entry : same_strand){
        if(lis(entry.second)>=lis_threshold){
            overlaps.emplace_back(entry.first,true);
        }
    }

    for(auto &entry : different_strand){
        if(lis(entry.second)>=lis_threshold){
            overlaps.emplace_back(make_pair(entry.first,false));
        }
    }
    return overlaps;
};

// Binary search (note boundaries in the caller)
int cell_index(std::vector<int> &v, int l, int r, int key) {
    while (r-l > 1) {
        int m = l + (r-l)/2;
        if (v[m] >= key)
            r = m;
        else
            l = m;
    }

    return r;
}

int lis(const std::vector<int> &v) {
    int v_size = v.size();
    if (v_size == 0)
        return 0;

    std::vector<int> tail(v_size, 0);
    int length = 1; // always points empty slot in tail

    tail[0] = v[0];
    for (size_t i = 1; i < v_size; i++) {
        if (v[i] < tail[0])
            // new smallest value
            tail[0] = v[i];
        else if (v[i] > tail[length-1])
            // v[i] extends largest subsequence
            tail[length++] = v[i];
        else
            // v[i] will become end candidate of an existing subsequence or
            // Throw away larger elements in all LIS, to make room for upcoming grater elements than v[i]
            // (and also, v[i] would have already appeared in one of LIS, identify the location and replace it)
            tail[cell_index(tail, -1, length-1, v[i])] = v[i];
    }

    return length;
}

