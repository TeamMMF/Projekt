//
// Created by filip on 07.12.17..
//
#include "ProjectConfig.h"
#include "Common.hpp"
#include <iostream>
#include "FASTASampleClass.cpp"
#include "bioparser/bioparser.hpp"

#define WINDOW_DEFAULT 5
#define KMER_DEFAULT 3


void show_usage(string arg) {
    fprintf(stdout, "%s Version %d.%d\n",
            arg.c_str(),
            SequenceOverlaping_VERSION_MAJOR,
            SequenceOverlaping_VERSION_MINOR
    );
    fprintf(stdout, "Usage: %s <option(s)> <inputfile>\n", arg.c_str());
    fprintf(stdout, "Available options are:\n");
    fprintf(stdout, "-h, --help\t\t\tshow usage instructions.\n");
    fprintf(stdout, "-w, --window\t\t\tspecify the window size, default is %d\n", WINDOW_DEFAULT);
    fprintf(stdout, "-k, --kmer\t\t\tspecify the k-mer size, default is %d\n", KMER_DEFAULT);

}

int main(int argc, char const *argv[]) {

    int k = KMER_DEFAULT;
    int w = WINDOW_DEFAULT;

    string sequence_file_path;

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
                if (i + 1 >= argc) { // Make sure we aren't at the end of argv!
                    std::cerr << "option --window takes one argument." << std::endl;
                    return 1;
                }
                w = atoi(argv[++i]);

            } else if ((arg == "-k") || (arg == "--kmer")) {
                if (i + 1 >= argc) { // Make sure we aren't at the end of argv!
                    std::cerr << "option --kmer takes one argument." << std::endl;
                    return 1;
                }
                k = atoi(argv[++i]);
            }else if (i + 1 == argc){
                sequence_file_path = argv[i];
            }else{
                show_usage(argv[0]);
                return 1;
            }
        }
    }

    if(sequence_file_path.empty()){
        show_usage(argv[0]);
    }


    vector<unique_ptr < FASTASampleClass>> fasta_reads;
    auto fasta_reader = bioparser::createReader<FASTASampleClass, bioparser::FastaReader>(sequence_file_path);
    fasta_reader->read_objects(fasta_reads, static_cast<uint64_t>(-1));

    vector<string> sequences;
    sequences.reserve(fasta_reads.size());

    for (auto &fasta_read : fasta_reads) {
        sequences.push_back(fasta_read->get_data());
    }
    cout << sequence_file_path<<endl;
    cout << w<<endl;
    cout << k<<endl;
    auto lookup = indexTable(sequences, w, k);


    return 0;
}

