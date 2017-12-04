//
// Created by matep on 03/12/2017.
//


#include <stdint-gcc.h>

void LCS(const char *s1, const char *s2, int s1_l, int s2_l);

void print_matrix(const char *s1, const char *s2, uint8_t *matrix, int r, int c);

bool check_k_substring_match(const char* s1, const char* s2, int i, int j, int k);

void LCS_k(const char* s1, const char* s2, int s1_l, int s2_l, int k);