//
// Created by matep on 03/12/2017.
//

#include <string>
#include <iostream>
#include <iomanip>
#include <cstring>
#include "Dynamic.h"

using namespace std;

void LCS( const char* s1, const char* s2, int s1_l, int s2_l){
    int r = s1_l + 1;
    int c = s2_l + 1;

    uint8_t *matrix = new uint8_t[r*c];

    memset(matrix, 0, r*c);

    print_matrix(matrix, r, c);

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
        print_matrix(matrix, r, c);
    }

    print_matrix(matrix, r, c);
    delete(matrix);
}

void print_matrix(uint8_t *matrix, int r, int c){
    for(int i = 0; i < r; i++){
        for(int j = 0; j < c; j++){
            cout << setw(4) << +matrix[i*r + j];
        }
        cout << endl;
    }
    cout << endl;
}