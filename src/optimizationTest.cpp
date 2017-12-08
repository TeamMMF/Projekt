//
// Created by matep on 08/12/2017.
//


#include "Common.hpp"

using namespace std;

int main(void){
    std::string s = "AACCTTGGAACCGGTTV";
    uint8_t k = 3;
    char** kmers;
    find_kmers(s.c_str(), 3, &kmers, s.length() - k + 1);
    printf("VRATIO SE IZ FUNKCIJE\n");
    for(int i = 0; i < s.length() - k + 1 ; i++){
        printf("%s %d\n", &(kmers[i]));
        free(kmers[i]);
    }

    //char* rev_comp = (char*) malloc(s.length() + 1);
    //find_reverse_complement(s.c_str(), s.length(), rev_comp);
    //int i = 0;
    //printf("%s", rev_comp);

    if(kmers != NULL){
        free(kmers);
    } else {
        printf("VEC JE OSLOBODEN");
    }
    //free(rev_comp);
    //uint16_t w = 5;
    //minimizer *minimizers = new minimizer[s.length() - k - w + 2];
    //find_minimizers3(s.c_str(),(uint16_t) s.length(), w, k, minimizers, s.length() - k - w + 2);
    //delete[] minimizers;
}