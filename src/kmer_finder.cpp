//
// Created by matep on 07/11/2017.
//

#include "ProjectConfig.h"

#include <string>
#include <vector>
#include <iostream>
#include <unordered_set>
#include <set>

#include <sstream>
#include <memory>


#include "FASTASampleClass.cpp"
#include "bioparser/bioparser.hpp"

using namespace std;

vector<string> find_kmer(int k, string s)
{
    vector<string> v;
    v.reserve(s.length() - k + 1);
    for(int i = 0, len = s.length() - k; i <= len; i++){
        v.push_back(s.substr(i, k));
    }
    return v;
}

size_t minimizer_hash(string s){
    return hash<string>{}(s);
}

vector<string> find_minimizers(int w, int k, string s)
{

    vector<string> kmers = find_kmer(k, s);
    vector<string> minimizers;
    minimizers.reserve(kmers.size() - w + 1);

    pair<size_t, string> * hashed_kmers = new pair<size_t, string>[kmers.size()];

    for(int i = 0; i < kmers.size(); i++){
        hashed_kmers[i] = make_pair(minimizer_hash(kmers[i]), kmers[i]);
    }

    for(int i = 0, len = kmers.size() - w + 1; i < len; i++){      //prvi prozor cijeli hashirat, sljedeci samo zadnju ktorku
        for(int j = 0; j < w; j++){
            string new_kmer = s.substr(i + j, k);
            size_t new_hash = minimizer_hash(new_kmer);
            set.insert(make_pair(new_hash, new_kmer));
        }
        min.push_back((*set.rend()).second);
    }
    return min;
}

int main(int argc, char const *argv[])
{

    string project_root(PROJECT_ROOT);

    vector<unique_ptr<FASTASampleClass>> fasta_objects;
    auto fasta_reader = bioparser::createReader<FASTASampleClass, bioparser::FastaReader>(project_root + "src/resources/sample.fasta");

    fasta_reader->read_objects(fasta_objects, static_cast<uint64_t>(-1));

    for (auto &fasta_object : fasta_objects) {

        //vector<string> vector = find_kmer(100, (*fasta_object).get_data());
        vector<string> vector = find_minimizers(20, 100, (*fasta_object).get_data());
        for(int i = 0, len = vector.size(); i < len; i++){
            cout << vector.at(i) << endl;
        }
        break;
    }

    /*
    vector<string> vector = find_kmer(3, s);
    for(int i = 0, len = vector.size(); i < len; i++)
    {
        cout << vector.at(i) << endl;
    }
     */

    return 0;
}

