#include <vector>
#include <string>

using namespace std;

typedef tuple<string, int, string> triplet;


char complement(char c);

std::vector<std::string> find_kmer(int k, std::string s);

size_t minimizer_hash(string s);

vector<triplet> find_minimizers(int w, int k, string s);