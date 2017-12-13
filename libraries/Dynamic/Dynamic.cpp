//
// Created by matep on 03/12/2017.
//
#include "Dynamic.h"
#include "Common.hpp"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <algorithm>
#include <unordered_map>
#include <climits>
#include <CustomTypes.h>


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
                     unordered_multimap<uint64_t, int>* seq2_hash_to_index,
                     minimizer* seq2_mins_sorted){
    vector<int> lis_arr;

    for(int i = 0; i < seq1_mins_size; i++){
        auto range = seq2_hash_to_index->equal_range(seq1_mins_sorted[i].hash);
        int min_diff = 0;
        int final = 0;
        for (auto it = range.first; it != range.second; ++it) {
            int cur_diff = it->second - seq1_mins_sorted->index;
            if(seq2_mins_sorted[it->second].rev == seq1_mins_sorted[i].rev
               && cur_diff < min_diff){
                final = it->second;
                min_diff = cur_diff;
            }
        }
        lis_arr.push_back(final);
    }

    return lis(&lis_arr[0], lis_arr.size());

}


int lis( int* a, int N ) {
  int *best, i, j, max = INT_MIN;

  best = (int*) malloc ( sizeof( int ) * N );
  for ( i = 0; i < N; i++ ) best[i] = 1;
 
  // best is an array initialized to all 1s becasue thats the shortest possible answer,
  // which would happen when all numbers in the array are always decreasing. We will use
  // this array to keep track of the best answer for up to each position as we traverse 
  // through the array below.

  // travese from 2nd element till the end
  // this is the right boundary for each iteration
    for ( i = 1; i < N; i++ )

      // traverse from 1st element until i'th element
      for ( j = 0; j < i; j++ )
      {
          // if:
          // a[j] is less than a[i] (number to the left of a[j])
          //     and
          // this increases the best count for that element
          if ( a[i] > a[j] && best[i] < best[j] + 1 )
          {
            // make a record of this
            best[i] = best[j] + 1;

            // if this best is greater than max, make a record of that
            if(max < best[i])
                  max = best[i];
           }         
       }
    // free
    free( best );

   // return
   return max;
}