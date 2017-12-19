//
// Created by filip on 19.12.17..
//

#include "ProjectConfig.h"
#include "Common.hpp"
#include "FASTARead.h"
#include "Dynamic.h"
#include <CustomTypes.h>
#include <chrono>
#include <lcskpp.h>
#include "bioparser/bioparser.hpp"
#include <algorithm>

#define WINDOW_DEFAULT 5
#define KMER_DEFAULT 15


const char *progress = "-\\|/";

bool lis_threshold(int result,int l1, int l2);


void show_usage(string arg) {
    fprintf(stdout, "%s Version %d.%d\n",
            arg.c_str(),
            SequenceOverlaping_VERSION_MAJOR,
            SequenceOverlaping_VERSION_MINOR
    );
    fprintf(stdout, "Usage: %s <option(s)> <readfile> <resultsfile>\n", arg.c_str());
    fprintf(stdout, "Available options are:\n");
    fprintf(stdout, "-h, --help\t\t\tshow usage instructions.\n");

}

void report_status(const char* operation, int curr, long total) {
    long ratio = 100*curr/total;
    fprintf(stdout,"\r%s [%c] %ld%c",operation, progress[curr%4],ratio, '%');
    fflush(stdout);
}

bool lis_threshold(int result, int l1, int l2) {
    return result > 9;
}

// Mislio si izvest veceg Matu Paulinovica... Take a seat Skywalker
int main(int argc, char const *argv[]) {

    uint32_t k = KMER_DEFAULT;
    uint32_t w = WINDOW_DEFAULT;

    string read_file_path;
    string result_file_path;

    if (argc < 3) {
        show_usage(argv[0]);
        return 1;
    }
    read_file_path = argv[1];
    result_file_path = argv[2];

    if (read_file_path.empty()) {
        show_usage(argv[0]);
        return 1;
    }
    // čita datoteku i sprema svako očitanje u poseban objekt razreda FASTARead
    vector<unique_ptr<FASTARead>> fasta_reads;
    auto fasta_reader1 = bioparser::createReader<FASTARead, bioparser::FastaReader>(read_file_path);
    fasta_reader1->read_objects(fasta_reads, static_cast<uint64_t>(-1));
    long number_of_reads = fasta_reads.size();
    printf("Reading file - Done\n");

    unordered_map<uint64_t, vector<hashMinPair2>> lookup_map; // hash minimizera -> minimizeri svih sekvenci poredani po indeksu uzlazno
    //std::vector<std::vector<uint64_t >> mins_in_order; // id sekvence -> poredani minimizeri sekvence po indeksu
    std::vector<std::vector<minimizer>> mins_in_order(number_of_reads);

    printf("Colecting data [-]");
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    for (int i=0; i<number_of_reads; i++){
        report_status("Collecting data",i, number_of_reads);
        process_sequence4(fasta_reads[i]->get_data(),
                          fasta_reads[i]->get_data_length(),
                          i,
                          w,
                          k,
                          mins_in_order);
    }

    fill_lookup_table(mins_in_order, lookup_map);
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    printf("\rCollecting data - Finished in %ld seconds\n",chrono::duration_cast<chrono::seconds>( t2 - t1 ).count());

    printf("Preparing data for processing.\n");
    fflush(stdout);
    chrono::high_resolution_clock::time_point t3 = chrono::high_resolution_clock::now();
    sort_by_indices(lookup_map);
    chrono::high_resolution_clock::time_point t4 = chrono::high_resolution_clock::now();
    printf("Data prepared in %ld seconds", chrono::duration_cast<chrono::seconds>( t4 - t3 ).count());
    fflush(stdout);
    printf("\nComparing sequences [-]");
    chrono::high_resolution_clock::time_point t5 = chrono::high_resolution_clock::now();
    FILE* output = fopen("out.paf","w");
    /*
    for (int i = 0; i < number_of_reads; ++i) {
        report_status("Comparing sequences",i, number_of_reads);
        vector<pair<int, bool>> result = find_overlaps_by_LIS(i,mins_in_order[i],lookup_map,6);
        for(auto res : result){
            fprintf(output, "%s\t%d\t%d\t%d\t%c\t%s\t%d\n",
                    fasta_reads[i] -> get_name(),
                    fasta_reads[i] -> get_data_length(),
                    0,
                    0,
                    res.second ? '+' : '-',
                    fasta_reads[res.first] -> get_name(),
                    fasta_reads[res.first] -> get_data_length()
            );
        }
    }
     */
    chrono::high_resolution_clock::time_point t6 = chrono::high_resolution_clock::now();
    printf("\rComparing sequences - Finished in %ld seconds.\n", chrono::duration_cast<chrono::seconds>( t6 - t5 ).count());
    printf("Total execution time: %ld seconds\n", chrono::duration_cast<chrono::seconds>( t6 - t1 ).count());

    return 0;
}