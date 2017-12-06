//
// Created by matep on 04/12/2017.
//

#include "Dynamic.h"
#include <iostream>
#include <unordered_map>
#include <functional>


int main(){
    LCS("AAA", "ABA", 3, 3);
    LCS("XMJYAUZ", "MZJAWXU", 7, 7);
    LCS_k("aaaaaaaa", "aaaaaaaa", 8, 8, 2);
    LCS_k("aabbccdd", "bbaaddcc", 8, 8, 2);
    LCS_k("ATTAT", "CTATAGAGTA", 5, 10, 2);
    LCS_k("ATTATG", "CTATAGAGTA", 6, 10, 2);

    generate_match_pairs("ATTAT", "CTATAGAGTA", 2);
    std::cout<< "new pairs" << std::endl;
    generate_match_pairs("ATTATG", "CTATAGAGTA", 2);
    std::cout<< "new pairs" << std::endl;
    generate_match_pairs("TATGATA", "TATGATA", 3);
    std::cout << unsigned(LCS_kpp("TATGATA", "TATGATA", 3)) << std::endl;
    std::cout << unsigned(LCS_kpp("ATCTAG", "ATCTAG", 3)) << std::endl;

    std::cout << unsigned(LCS_kpp("ATCTTTTAG", "ATCGGGTAG", 3)) << std::endl;

    check_precedence(std::make_tuple(1,1,true), std::make_tuple(2,2,true));
    /*
    std::unordered_map<matchPair, int, std::function<size_t(std::tuple<int,int,bool>)>> ordering(10, tuple_hash);
    std::cout << tuple_hash(std::make_tuple(1,1,true)) << std::endl;
    ordering.emplace(std::make_tuple(1,1,true), 2);
    std::cout << ordering.find(std::make_tuple(1,1,true))->second;
    */
    return 0;
}