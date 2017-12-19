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


int lds(vector<int> vector);

using namespace std;


typedef struct{
    uint64_t  seq_id;
    int index;
    bool reverse;
} tir;

void LCS( const char* s1, const char* s2, int s1_l, int s2_l){
    int r = s1_l + 1;
    int c = s2_l + 1;

    uint8_t *matrix = new uint8_t[r*c];

    memset(matrix, 0, r*c);

    //print_matrix(matrix, r, c);

    for(int i = 1; i < r; i++){
        for(int j = 1; j < c; j++){
            if(s1[i-1] == s2[j-1]){
                matrix[i*r + j] = matrix[(i-1)*r + (j-1)] + 1;
                continue;
            }

            uint8_t left = matrix[(i-1)*r + j];
            uint8_t up = matrix[i*r + (j-1)];
            if(left < up){
                matrix[i*r + j] = up;
            }
            else{
                matrix[i*r + j] = left;
            }
        }
        //print_matrix(matrix, r, c);
    }

    print_matrix(s1, s2, matrix, r, c);

    string s = traceback_LCS(matrix, r, c, s1);
    cout << "LCS: " << s << endl;
    delete[] matrix;
}

string traceback_LCS(uint8_t *matrix, int r, int c, const char* s1){
    std::string lcs;
    int i = r - 1;
    int j = c - 1;
    while(j > 0 || i > 0){
        if(j == 0){
            i--;
            continue;
        }
        else if(i == 0){
            j--;
            continue;
        }

        uint8_t up = matrix[(i-1)*r + j];
        uint8_t left = matrix[i*r + (j-1)];
        if(up > left){
            --i;
        }
        else if (left > up){
            --j;
        }
        else{
            uint8_t diag = matrix[(i-1)*r + (j-1)];
            if(diag < matrix[i*r + j]){
                i--;
                j--;
                lcs += s1[i];
            }
            else {

            }
        }
    }
    std::reverse(lcs.begin(), lcs.end());
    return lcs;
}

void LCS_k(const char* s1, const char* s2, int s1_l, int s2_l, int k){
    int r = s1_l + 1;
    int c = s2_l + 1;

    uint8_t *matrix = new uint8_t[r*c];

    memset(matrix, 0, r*c);

    //print_matrix(matrix, r, c);

    for(int i = 1; i < r; i++){
        for(int j = 1; j < c; j++){
            if(check_k_substring_match(s1, s2, i, j, k)){
                matrix[i*c + j] = matrix[(i-k)*c + (j-k)] + k;
                cout << "MATCH PAIR (" << i << ", " << j << ")" << endl;
                continue;
            }

            uint8_t left = matrix[i*c + (j-1)];
            uint8_t up = matrix[(i-1)*c + j];
            if(left < up){
                matrix[i*c + j] = up;
            }
            else{
                matrix[i*c+ j] = left;
            }
        }
        //print_matrix(s1, s2, matrix, r, c);
    }

    print_matrix(s1, s2, matrix, r, c);
    delete[] matrix;
}

bool check_k_substring_match(const char* s1, const char* s2, int i, int j, int k){
    if(i < k || j < k) {
        return false;
    }

    for(int l = 0; l < k; l++){
        if(s1[i - l - 1] != s2[j - l - 1]){
            return false;
        }
    }

    return true;
}

void print_matrix(const char *s1, const char *s2, uint8_t *matrix, int r, int c){
    cout << "     ";
    for(int i = 0; i < c - 1; i++){
        cout << setw(4) << s2[i];
    }
    cout << endl;

    for(int i = 0; i < r; i++){
        if(i > 0){
            cout << s1[i - 1];
        }
        else {
            cout << " ";
        }
        for(int j = 0; j < c; j++){
            cout << setw(4) << +matrix[i*c + j];
        }
        cout << endl;
    }
    cout << endl;
}

vector<matchPair> generate_match_pairs(string s1, string s2, int k){

    unordered_multimap<uint64_t, int> kmer_hash;
    for(int i = 0, len = s2.length() - k ; i <= len; i++){
        kmer_hash.emplace(minimizer_hash(s2.substr(i, k)), i);

    }

    vector<matchPair> match_points;


    for(int i = 0, len = s1.length() - k; i <= len; i++){
        string sub = s1.substr(i, k);
        auto column_hits = kmer_hash.equal_range(minimizer_hash(sub));

        for(auto it = column_hits.first; it != column_hits.second; ++it){
            match_points.emplace_back(i, it->second , true);
            match_points.emplace_back(i + k, it->second + k, false);
        }

    }

    sort(match_points.begin(), match_points.end(), matchPair_comparator);
    //for(auto mp : match_points){
    //    cout << "(" << get<0>(mp) << ", " << get<1>(mp) << ", start = "<< get<2>(mp) << ")" << endl;
    //}

    return match_points;
}

bool matchPair_comparator(const matchPair a, const matchPair b){
    if(get<0>(a) < get<0>(b)) return true;
    if(get<0>(a) > get<0>(b)) return false;

    if(get<1>(a) < get<1>(b)) return true;
    if(get<1>(a) > get<1>(b)) return false;

    if( get<2>(a) &&  get<2>(b)) return true;
    if( get<2>(a) && !get<2>(b)) return false;
    if(!get<2>(a) &&  get<2>(b)) return true;
    if(!get<2>(a) && !get<2>(b)) return true;

    return true;
}

uint8_t max_vector(vector<uint8_t> vector){
    uint8_t max = vector[0];
    for(int i = 1, len = vector.size(); i < len; i++){
        uint8_t tmp = vector[i];
        if(tmp > max){
            max = tmp;
        }
    }

    return max;
}


size_t tuple_hash(std::tuple<int, int, bool> x){
    std::hash<int> int_hash;
    std::hash<bool> bool_hash;
    return (int_hash(std::get<0>(x)) ^ int_hash(std::get<1>(x)) + 13 * bool_hash(std::get<2>(x)));
}


uint64_t LCS_kpp(string s1, string s2, int k){
    vector<uint64_t> max_col_dp(s2.size() + 1);

    vector<matchPair> match_pairs = generate_match_pairs(s1, s2, k);

    unordered_map<matchPair, int, function<size_t(matchPair)>>  ordering(match_pairs.size() / 2, tuple_hash); // codes matchPair to index in vector
    unordered_map<int, int> row_to_vector_position(match_pairs.size()/2); // binds row in "matrix" to first matchPair position in match_pairs vector

    vector<uint64_t> dp(match_pairs.size()/2);

    int start_counter = 0;   //counts number of start events
    int n = 0;              //keeps track of position within match_pairs vector
    for(matchPair mp : match_pairs){

        if(get<2>(mp)) {             //match pair is start
            ordering.emplace(mp, start_counter);
            dp[start_counter] = k + max_between_indexes(max_col_dp, 0, get<1>(mp));     //line 6 in pseudocode
            start_counter++;

            auto t = row_to_vector_position.find(get<0>(mp));   //adds another entry to row_to_vector_position if one does not exist already
            if(t == row_to_vector_position.end()){
                row_to_vector_position.emplace(get<0>(mp), n);
            }
        }


        else {
            int n_pred = -1;
            int n_mp = ordering.find(make_tuple(get<0>(mp) - k, get<1>(mp) - k, true))->second;  //find the index of the current start matchPoint in the dp vector
            auto m_it = row_to_vector_position.find(get<0>(mp) - k - 1);
            if(m_it != row_to_vector_position.end()) {                      //if there exists a start in that row

                int m = m_it->second;

                for (; m != n; m++) {
                    matchPair pred_candidate = match_pairs[m];

                    bool precedence = false;
                    if (!get<2>(pred_candidate)) {                              //if event is not start continue
                        continue;
                    }
                    else if(get<0>(pred_candidate) > (get<0>(mp) - k + 1)){     //if event has gone to the next row
                        break;
                    }

                    int p_diff = (get<0>(mp) - k) - (get<1>(mp) - k);               //calculates if mp continues the other event
                    int g_diff = get<0>(pred_candidate) - get<1>(pred_candidate);
                    int i_diff = (get<0>(mp) - k) - get<0>(pred_candidate);

                    precedence = (p_diff == g_diff) && (i_diff == 1);

                    if (precedence) {                                               //if the event continues
                        /*if((get<0>(mp) - get<0>(pred_candidate) == k) && (get<1>(mp) - get<1>(pred_candidate) == k)){
                            break;
                        }*/
                        n_pred = ordering.find(pred_candidate)->second;      //find the index of the matchPoint which mp continues
                        break;
                    }

                }
            }

            if(n_pred != -1){
                dp[n_mp] = max(dp[n_mp], dp[n_pred] + 1);
            }

            max_col_dp[get<1>(mp)] = max(max_col_dp[get<1>(mp)], dp[n_mp]);
        }

        n++;
    }

    return max_between_indexes(dp, 0, dp.size() - 1);
}
bool check_precedence(matchPair p, matchPair g){
    if(!get<2>(p)){
        return false;
    }

    int p_diff = get<0>(p) - get<1>(p);
    int g_diff = get<0>(g) - get<1>(g);
    int i_diff = get<0>(p) - get<0>(g);

    return (p_diff == g_diff) && (i_diff == 1);
}

uint64_t max_between_indexes(vector<uint64_t> array, int start, int end){ // both inclusive
    uint64_t max = 0;
    uint64_t test;

    for(int i = start; i <= end; i++){
        test = array[i];
        if(max < test){
            max = test;
        }
    }

    return max;
}


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
                                                  vector<uint64_t>& minimizer_hashes,
                                                   unordered_map<uint64_t, vector<hashMinPair2>>&  minimizers_for_hash,
                                                  int lis_threshold){
    unordered_map<uint64_t, vector<int>> same_strand;
    unordered_map<uint64_t, vector<int>> different_strand;

    for(auto h : minimizer_hashes){
        auto matches = minimizers_for_hash.find(h);
        if(matches == minimizers_for_hash.end())
            continue;
        for(auto match : matches->second){
            if(match.seq_id <= query_id){
                continue;
            }
            if(match.rev){
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

