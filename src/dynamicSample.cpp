//
// Created by matep on 04/12/2017.
//

#include "Dynamic.h"

int main(){
    LCS("AAA", "ABA", 3, 3);
    LCS("XMJYAUZ", "MZJAWXU", 7, 7);
    LCS_k("aaaaaaaa", "aaaaaaaa", 8, 8, 2);
    LCS_k("aabbccdd", "bbaaddcc", 8, 8, 2);
    LCS_k("ATTAT", "CTATAGAGTA", 5, 10, 2);
    LCS_k("ATTATG", "CTATAGAGTA", 6, 10, 2);
    uint64_t  d = test();
//    generate_match_pairs("aabbccdd", "bbaaddcc", 2);
    return 0;
}