//
// Created by matep on 03/12/2017.
//

#include <string>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <algorithm>
#include "Common.hpp"
#include "Dynamic.h"

using namespace std;

void LCS( const char* s1, const char* s2, int s1_l, int s2_l){
    int r = s1_l + 1;
    int c = s2_l + 1;

    uint8_t *matrix = new uint8_t[r*c];

    memset(matrix, 0, r*c);
    char a = complement('A');

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
    delete(matrix);
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
                matrix[i*r + j] = matrix[(i-k)*r + (j-k)] + k;
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
    delete(matrix);
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
    for(int i = 0; i < r - 1; i++){
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
            cout << setw(4) << +matrix[i*r + j];
        }
        cout << endl;
    }
    cout << endl;
}