#include <vector>
#include <string>

using namespace std;

typedef tuple<string, int, string> triplet;

typedef tuple<uint64_t, int, int> triple;

char complement(char c);

std::vector<std::string> find_kmer(int k, std::string s);

uint64_t minimizer_hash(string s);

uint64_t invertible_minimizer_hash(uint64_t x);

uint64_t invertible_minimizer_hash_inverse(uint64_t x);

vector<triplet> find_minimizers(int w, int k, string s);

vector<tuple<uint64_t, int, int>> find_minimizers2(int w, int k, string s);

size_t find_hash_value(char c);