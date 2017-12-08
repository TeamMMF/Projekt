#include <vector>
#include <string>
#include <functional>
#include <unordered_map>
#include <CustomTypes.h>

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

std::unordered_multimap<uint64_t, tuple<string, int, int>, function<size_t( uint64_t)>> indexSequences(vector<string> sequences, int w, int k);
std::unordered_multimap<uint64_t, hashEntry, function<size_t( uint64_t)>> indexSequences(string sequences, int w, int k);

const int hash_width = 20;
const int offset_width = 6;

typedef tuple<uint64_t, tuple<string, int, int>> hashMinPair;

std::vector<hashMinPair> indexTable(vector<string> sequences, int w, int k);

typedef tuple<string, int, int, int> minimizer_hit;

bool hit_comparator(const minimizer_hit a,const minimizer_hit b);

typedef struct {
    int target_min_index;
    int target_max_index;
    int query_min_index;
    int query_max_index;
    //string query_name;
    bool reverse;
} mapInfo;

vector<mapInfo> map_minimizers(unordered_multimap<uint64_t, tuple<string, int, int>, function<size_t(uint64_t)>> lookup_table, string query_sequence, int w, int k, int epsilon);
