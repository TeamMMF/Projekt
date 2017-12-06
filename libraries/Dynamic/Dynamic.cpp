//
// Created by matep on 03/12/2017.
//

#include <string>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <algorithm>
#include <vector>
#include <unordered_map>
#include "Dynamic.h"
#include "Common.hpp"

using namespace std;

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

uint64_t computeJustMin(string s){
    return minimizer_hash(s);
}

void generate_match_pairs(string s1, string s2, int k){

    unordered_multimap<uint64_t, int> kmer_hash;
    for(int i = 0, len = s2.length() - k ; i <= len; i++){
        kmer_hash.emplace(minimizer_hash(s2.substr(i, k)), i);

    }

    vector<matchPair> match_points;


    for(int i = 0, len = s1.length() - k; i <= len; i++){
        string sub = s1.substr(i, k);
        auto column_hits = kmer_hash.equal_range(minimizer_hash(sub));

        for(auto it = column_hits.first; it != column_hits.second; ++it){
            match_points.emplace_back(i + 1, it->second + 1, true);
            match_points.emplace_back(i+k, it->second + k, false);
        }

    }


    sort(match_points.begin(), match_points.end(), matchPair_comparator);
    for(auto mp : match_points){
        cout << "(" << get<0>(mp) << ", " << get<1>(mp) << ", start = "<< get<2>(mp) << ")" << endl;
    }

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
