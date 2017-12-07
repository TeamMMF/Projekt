//
// Created by filip on 07.12.17..
//
#include "ProjectConfig.h"
#include "Common.hpp"
#include "Dynamic.h"
#include <iostream>
#include <cstdio>
#include "FASTASampleClass.cpp"
#include "bioparser/bioparser.hpp"

#define WINDOW_DEFAULT 5
#define KMER_DEFAULT 3
#define EPSILON_DEFAULT 1


void show_usage(string arg) {
    fprintf(stdout, "%s Version %d.%d\n",
            arg.c_str(),
            SequenceOverlaping_VERSION_MAJOR,
            SequenceOverlaping_VERSION_MINOR
    );
    fprintf(stdout, "Usage: %s <option(s)> <readfile> <referencefile>\n", arg.c_str());
    fprintf(stdout, "Available options are:\n");
    fprintf(stdout, "-h, --help\t\t\tshow usage instructions.\n");
    fprintf(stdout, "-w, --window\t\t\tspecify the window size, default is %d\n", WINDOW_DEFAULT);
    fprintf(stdout, "-k, --kmer\t\t\tspecify the k-mer size, default is %d\n", KMER_DEFAULT);
    fprintf(stdout, "-e, --kmer\t\t\tspecify the epsilon value, default is %d\n", EPSILON_DEFAULT);

}

int main(int argc, char const *argv[]) {

    int k = KMER_DEFAULT;
    int w = WINDOW_DEFAULT;
    int eps = EPSILON_DEFAULT;

    string read_file_path;
    string reference_file_path;

    if (argc < 2) {
        show_usage(argv[0]);
        return 1;
    }
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-h") || (arg == "--help")) {
            show_usage(argv[0]);
            return 0;
        } else {

            if ((arg == "-w") || (arg == "--window")) {
                if (i + 1 >= argc) {
                    std::cerr << "option --window takes one argument." << std::endl;
                    return 1;
                }
                w = atoi(argv[++i]);

            } else if ((arg == "-k") || (arg == "--kmer")) {
                if (i + 1 >= argc) {
                    std::cerr << "option --kmer takes one argument." << std::endl;
                    return 1;
                }
                k = atoi(argv[++i]);

            } else if ((arg == "-e") || (arg == "--epsilon")) {
                if (i + 1 >= argc) {
                    std::cerr << "option --epsilon takes one argument." << std::endl;
                    return 1;
                }
                eps = atoi(argv[++i]);

            }else if (i + 2 == argc){
                read_file_path = argv[i];
            }else if (i + 1 == argc){
                reference_file_path = argv[i];
            }else{
                show_usage(argv[0]);
                return 1;
            }
        }
    }

    if(read_file_path.empty()){
        show_usage(argv[0]);
        return 1;
    }

    vector<unique_ptr < FASTASampleClass>> fasta_reads;
    auto fasta_reader1 = bioparser::createReader<FASTASampleClass, bioparser::FastaReader>(read_file_path);
    fasta_reader1->read_objects(fasta_reads, static_cast<uint64_t>(-1));

    vector<string> read_sequences;
    read_sequences.reserve(fasta_reads.size());

    for (auto &fasta_read : fasta_reads) {
        read_sequences.push_back(fasta_read->get_data());
    }

    vector<unique_ptr < FASTASampleClass>> fasta_reference;
    auto fasta_reader2 = bioparser::createReader<FASTASampleClass, bioparser::FastaReader>(reference_file_path);
    fasta_reader2->read_objects(fasta_reference, static_cast<uint64_t>(-1));

    vector<string> target_sequences;
    read_sequences.reserve(fasta_reference.size());
    vector<string> box;
    box.reserve(fasta_reference.size());

    for (auto &fasta_ref : fasta_reference) {
        box.push_back(fasta_ref->get_data());
        target_sequences.push_back(fasta_ref->get_data());
    }

    cout << read_file_path<<endl;
    cout << reference_file_path<<endl;
    cout << w<<endl;
    cout << k<<endl;
    for(auto &seq : read_sequences){
        vector<mapInfo> map_info = map_minimizers(indexSequence(box,w,k), seq, k, w, eps);
        int len = seq.size();
        for(auto &mapinf : map_info){
            uint64_t similarity = LCS_kpp(seq.substr(mapinf.target_min_index-len/2,
                                                     mapinf.target_max_index-mapinf.target_min_index+1),
                                          seq,
                                          k);
            fprintf(stdout,"%d\t%d\t%d\t%d\t%d\t\t%ld",
                    mapinf.query_min_index,
                    mapinf.query_max_index,
                    mapinf.target_min_index,
                    mapinf.target_max_index,
                    mapinf.reverse,
                    similarity);
        }
    }

    return 0;
}

