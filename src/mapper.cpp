//
// Created by filip on 07.12.17..
//
#include "ProjectConfig.h"
#include "Common.hpp"
#include "FASTARead.h"
#include "Dynamic.h"
#include <iostream>
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
    fprintf(stdout, "Usage: %s <option(s)> <readfile>\n", arg.c_str());
    fprintf(stdout, "Available options are:\n");
    fprintf(stdout, "-h, --help\t\t\tshow usage instructions.\n");
    fprintf(stdout, "-w, --window\t\t\tspecify the window size, default is %d\n", WINDOW_DEFAULT);
    fprintf(stdout, "-k, --kmer\t\t\tspecify the k-mer size, default is %d\n", KMER_DEFAULT);
   // fprintf(stdout, "-e, --kmer\t\t\tspecify the epsilon value, default is %d\n", EPSILON_DEFAULT);

}

int main(int argc, char const *argv[]) {

    uint32_t k = KMER_DEFAULT;
    uint32_t w = WINDOW_DEFAULT;
    uint32_t eps = EPSILON_DEFAULT;

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

            } else if (i + 1 == argc) {
                read_file_path = argv[i];
            } else {
                show_usage(argv[0]);
                return 1;
            }
        }
    }

    if (read_file_path.empty()) {
        show_usage(argv[0]);
        return 1;
    }

    // čita datoteku i sprema svako očitanje u poseban objekt razreda FASTARead
    vector<unique_ptr<FASTARead>> fasta_reads;
    auto fasta_reader1 = bioparser::createReader<FASTARead, bioparser::FastaReader>(read_file_path);
    fasta_reader1->read_objects(fasta_reads, static_cast<uint64_t>(-1));
    long number_of_reads = fasta_reads.size();


    // polje koje mapira indeks sekvence -> mapa minimizera
    auto min_hash_to_index = (unordered_multimap<uint64_t, int>*)
                    malloc(sizeof(unordered_multimap<uint64_t, int>)*fasta_reads.size());
    // polje koje mapira indeks sekvence -> adresa poredanog polja minimizera
    auto mins_in_order = (minimizer**) malloc(number_of_reads*sizeof(minimizer**));
    //polje koje mapira indeks sekvence -> velicina poredanog polja minimizera
    auto mins_number = (uint32_t*) malloc(number_of_reads*sizeof(uint32_t));

    // napuni inicijalizirana polja tako da indeks pojedinog ocitanja (iz fasta_reads) odgovara indeksu njegove mape minimizera
    // i indeksu polja njegovih poredanih minimizera
    for (int i=0; i<number_of_reads; i++){
        process_sequence(fasta_reads[i]->get_data(),
                         fasta_reads[i]->get_data_length(),
                         w,
                         k,
                         min_hash_to_index+i,
                         mins_in_order+i,
                         mins_number+i);
    }

    for(int i = 0; i < number_of_reads; i++){
        for (int j = i+1; j < number_of_reads; ++j) {
            int lis_result = compare_with_lis(*(mins_in_order+i),
                                              mins_number[i],
                                              min_hash_to_index+j,
                                              *(mins_in_order+j));
            //todo odluciti sto cemo dalje
        }
    }

    return 0;
}
