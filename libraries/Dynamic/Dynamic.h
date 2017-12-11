//
// Created by matep on 03/12/2017.
//
#ifndef Dynamic
#define Dynamic
#include <tuple>
#include <string>
#include <vector>
#include <cstdint>

void LCS(const char *s1, const char *s2, int s1_l, int s2_l);

void print_matrix(const char *s1, const char *s2, uint8_t *matrix, int r, int c);

bool check_k_substring_match(const char* s1, const char* s2, int i, int j, int k);

void LCS_k(const char* s1, const char* s2, int s1_l, int s2_l, int k);

std::string traceback_LCS(uint8_t *matrix, int r, int c, const char* s1);

typedef std::tuple<int, int, bool> matchPair;

std::vector<matchPair> generate_match_pairs(std::string s1, std::string s2, int k);

bool matchPair_comparator(const matchPair a, const matchPair b);

uint64_t max_between_indexes(std::vector<uint64_t> array, int start, int end);

bool check_precedence(matchPair p, matchPair g);

uint64_t LCS_kpp(std::string s1, std::string s2, int k);

size_t tuple_hash(std::tuple<int, int, bool> x);

int lis(int *a, int N);

#endif