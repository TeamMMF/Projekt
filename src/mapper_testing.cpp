//
// Created by filip on 07.12.17..
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

/*
const char *progress = "-\\|/";

bool lis_threshold(int result,int l1, int l2);

void report_status(const char string[16], int i, long reads);

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
    fprintf(stdout,"Reading file - Done\n");

    // polje koje mapira indeks sekvence -> mapa minimizera;
    unordered_map<int,unordered_multimap<uint64_t, int>> min_hash_to_index(number_of_reads);
    // polje koje mapira indeks sekvence -> adresa poredanog polja minimizera
    auto mins_in_order = (minimizer**) malloc(number_of_reads*sizeof(minimizer*));
    //polje koje mapira indeks sekvence -> velicina poredanog polja minimizera
    auto mins_number = (uint32_t*) malloc(number_of_reads*sizeof(uint32_t));

    // napuni inicijalizirana polja tako da indeks pojedinog ocitanja (iz fasta_reads) odgovara indeksu njegove mape minimizera
    // i indeksu polja njegovih poredanih minimizera
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    printf("Colecting data [-]");
    for (int i=0; i<number_of_reads; i++){
        report_status("Collecting data",i, number_of_reads);
        unordered_multimap<uint64_t, int> map(fasta_reads[i]->get_data_length() - w - k + 2);
        process_sequence(fasta_reads[i]->get_data(),
                         fasta_reads[i]->get_data_length(),
                         w,
                         k,
                         map,
                         mins_in_order+i,
                         mins_number+i);
        min_hash_to_index.emplace(i,map);
    }
    fprintf(stdout,"\rCollecting data - Done    ");
    fprintf(stdout,"\nComparing sequences [-]");
    FILE* output = fopen("out.paf","w");
    for(int i = 0; i < number_of_reads; i++){
        report_status("Comparing sequences",i, number_of_reads);
        for (int j = i+1; j < number_of_reads; ++j) {
            //printf(" %d,%d", i , j);
            pair<int,char> lis_result = compare_with_lis(mins_in_order[i],
                                              mins_number[i],
                                              min_hash_to_index.find(j)->second,
                                              mins_in_order[j]);
            if(!lis_threshold(lis_result.first,mins_number[i],
                                          mins_number[j])){
                continue;
            }
            //printf("END\n");
            fprintf(output, "%s\t%d\t%d\t%d\t%c\t%s\t%d\n",
                    fasta_reads[i] -> get_name(),
                    fasta_reads[i] -> get_data_length(),
                    0,
                    0,
                    lis_result.second,
                    fasta_reads[j] -> get_name(),
                    fasta_reads[j] ->get_data_length()
            );
        }
    }
    fprintf(stdout,"\rComparing sequences - Done     \n");

/* ZOVI ME MATE PAULINOVIC  VALGRIND

    fprintf(stdout,"\nComparing sequences naive [-]");
    FILE* output_n = fopen("out_naive.paf","w");
    for(int i = 0; i < number_of_reads; i++){
        report_status("Comparing sequences naive",i, number_of_reads);
        for (int j = i+1; j < number_of_reads; ++j) {
            int lis_same_strand = compare_with_lis(mins_in_order[i],
                                                   mins_number[i],
                                                   min_hash_to_index.find(j)->second,
                                                   mins_in_order[j],true);
            int lis_diff_strand = compare_with_lis(mins_in_order[i],
                                                   mins_number[i],
                                                   min_hash_to_index.find(j)->second,
                                                   mins_in_order[j],false);
            int lis_final;
            char strand;
            if(lis_diff_strand > lis_same_strand){
                lis_final = lis_diff_strand;
                strand ='-';
            } else {
                lis_final = lis_same_strand;
                strand = '+';
            }
            if(!lis_threshold(lis_final,mins_number[i],
                              mins_number[j])){
                continue;
            }
            fprintf(output_n, "%s\t%d\t%d\t%d\t%c\t%s\t%d\n",
                    fasta_reads[i] -> get_name(),
                    fasta_reads[i] -> get_data_length(),
                    0,
                    0,
                    strand,
                    fasta_reads[j] -> get_name(),
                    fasta_reads[j] ->get_data_length()
            );
        }
    }
    fprintf(stdout,"\rComparing sequences - Done     \n");*//*


    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    long duration = chrono::duration_cast<chrono::seconds>( t2 - t1 ).count();
    fprintf(stdout,"Completed in %ld seconds, results can be foud in the file %s.\n",duration,result_file_path.c_str());

    return 0;
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
    fprintf(stdout,"Reading file - Done\n");


    std::vector<std::vector<uint64_t >> mins_in_order;
    std::unordered_map<uint64_t , std::vector<hashMinPair2>> minimizer_hits;

    // napuni inicijalizirana polja tako da indeks pojedinog ocitanja (iz fasta_reads) odgovara indeksu njegove mape minimizera
    // i indeksu polja njegovih poredanih minimizera
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    printf("Colecting data [-]");
    printf("\nSizeof minimizer = %d\n", sizeof(minimizer));
    printf("Number of sequences to minimize: %ld\n", number_of_reads);

    uint64_t size = 0;
    uint64_t min_found = 0;
    for (int i=0; i<number_of_reads; i++){
        //report_status("Collecting data",i, number_of_reads);

        process_sequence3(fasta_reads[i]->get_data(),
                          fasta_reads[i]->get_data_length(),
                          i,
                          w,
                          k,
                          mins_in_order,
                          minimizer_hits);
        printf("READ %5d, minimizers found:%5d, map_size: %5d, seq_l: %5d\n",i, mins_in_order[i].size(), minimizer_hits.size(), fasta_reads[i]->get_data_length());
        size += mins_in_order[i].size() * sizeof(minimizer);
        min_found += mins_in_order[i].size();
        //if(i % 100 == 0){
        //    printf("USED BYTES SO FAR: %20llu\n",size);
        //}
    }
    printf("MINIMIZERS FOUND: %llu, TOTAL BYTE USAGE:%llu\n", min_found, size);
    fprintf(stdout,"\rCollecting data - Done    ");



    /*sort(minimizer_hits.begin(), minimizer_hits.end(), hashMinPair_comparator);
     *
    unordered_map<uint64_t, uint64_t> lookup_table;
    fill_lookup_table(&minimizer_hits, &lookup_table);

    printf("%d %d", minimizer_hits.size(), lookup_table.size());
    for(int i = 0; i < 100; i++){
        printf("%llu -> %llu\n", minimizer_hits[i].hash, lookup_table.find(minimizer_hits[i].hash)->second);
    }
    for(int i = 0; i < 100; i++){
        hashMinPair tmp = minimizer_hits[i];
        fprintf(stdout, "%3d, %20llu, %10u, %10u, %s\n", i, tmp.hash, tmp.seq_id, tmp.index, tmp.rev ? "True" : "False");
    }


    return 0;
}
 */