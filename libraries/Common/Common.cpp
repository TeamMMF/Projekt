#include <string>
#include <vector>
#include <tuple>
#include "Common.hpp"


using namespace std;

char complement(char c) {
    switch (c) {
        case 'A':
            return 'T';
        case 'T':
            return 'A';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        default:
            throw runtime_error("Invalid value.");
    }
}

vector<string> find_kmer(int k, string s) {
    vector<string> v;
    v.reserve(s.length() - k + 1);
    for (int i = 0, len = s.length() - k; i <= len; i++) {
        v.push_back(s.substr(i, k));
    }
    return v;
}

size_t minimizer_hash(string s) {
    return atoi(s.c_str());
}


std::vector<triplet> find_minimizers(int w, int k, string s) {

    vector<string> kmers = find_kmer(k, s);
    vector<triplet> minimizers;
    const unsigned long size = kmers.size();
    minimizers.reserve(size - w + 1);

    auto *hashed_kmers = new pair<size_t, string>[size];

    for (int i = 0; i < size; i++) {
        hashed_kmers[i] = make_pair(minimizer_hash(kmers[i]), kmers[i]);
    }

    for (int i = 0, len = size - w + 1; i < len; i++) {      //prvi prozor cijeli hashirat, sljedeci samo zadnju ktorku
        auto min = hashed_kmers[i];
        int minJ = i;
        for (int j = 1; j < w; j++) {
            auto current = hashed_kmers[i + j];
            if (min.first <= current.first) {
                continue;
            }
            min = current;
            minJ = i + j;
        }
        minimizers.emplace_back(make_tuple("asfda", minJ, min.second));
    }
    return minimizers;
}