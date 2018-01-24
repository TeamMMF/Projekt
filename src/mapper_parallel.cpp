//
// Created by filip on 19.12.17..
//

#include "ProjectConfig.h"
#include "Common.hpp"
#include "FASTARead.h"
#include "Dynamic.h"
#include <CustomTypes.h>
#include <chrono>
#include "bioparser/bioparser.hpp"
#include <algorithm>
#include <include/thread_pool/thread_pool.hpp>

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

void reduce_minimizers(uint32_t  number_of_minimizers,
                       std::unordered_map<uint64_t ,uint32_t>& min_occurences,
                       vector<uint64_t>& nogos){
    std::vector<std::pair<uint64_t, uint32_t>> min_occur(min_occurences.begin(), min_occurences.end());
    sort(min_occur.begin(), min_occur.end(), occurences_comparator);
    double acc = 0;
    for(int i = 0, len = min_occur.size(); i < len; i++){
        acc += min_occur[i].second / (double) number_of_minimizers;

        if(min_occur[i].second <= 34){
            break;
        }
        nogos.emplace_back(min_occur[i].first);
    }
}

uint32_t add_to_lookup_table(uint32_t seq_id,
                         std::vector<minim> &minimizers,
                         std::unordered_map<uint64_t,
                                 vector<hashMinPair3>>& map,
                         std::unordered_map<uint64_t ,uint32_t>& min_occurences){
    uint32_t  len = minimizers.size();
    for(uint32_t j = 0; j < len; j++){
        minim min = minimizers[j];
        auto it = map.find(min.hash);
        if (it == map.end()) {
            std::vector<hashMinPair3> vec;
            vec.emplace_back((hashMinPair3) {seq_id, min.index});
            map.emplace(min.hash, vec);
        } else {
            it->second.emplace_back((hashMinPair3) {seq_id, min.index});
        }

        min_occurences[min.hash]++;
    }
    return len;
}

void report_status(const char* operation, int curr, long total) {
    long ratio = 100*curr/total;
    fprintf(stdout,"\r%s [%c] %ld%c",operation, progress[curr%4],ratio, '%');
    fflush(stdout);
}

bool lis_threshold(int result, int l1, int l2) {
    return result > 9;
}

void lis_overlap_parallelization(int  query_id,
                                 vector<minim>& minimizer_hashes,
                                 unordered_map<uint64_t, vector<hashMinPair3>>&  lookup_map,
                                 int lis_threshold,
                                 vector<unique_ptr<FASTARead>>& fasta_reads,
                                 FILE* output,
                                 unordered_map<uint64_t,uint32_t >& occurences){

    vector<pair<int, bool>> result = find_overlaps_by_LIS_parallel(query_id,minimizer_hashes,lookup_map,lis_threshold,occurences);
    for(auto res : result){
        fprintf(output, "%s\t%d\t%d\t%d\t%c\t%s\t%d\n",
                fasta_reads[query_id] -> get_name(),
                fasta_reads[query_id] -> get_data_length(),
                0,
                0,
                res.second ? '+' : '-',
                fasta_reads[res.first] -> get_name(),
                fasta_reads[res.first] -> get_data_length()
        );
    }
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

    int max_concurrent_threads = argc > 3 ? stoi(argv[3]) : -1;
    if(max_concurrent_threads < 0){
        printf("Running on an unlimited number of concurrent threads.\n");
    }else{
        printf("Running on %d concurrent threads.\n",max_concurrent_threads);
    }

    // čita datoteku i sprema svako očitanje u poseban objekt razreda FASTARead
    vector<unique_ptr<FASTARead>> fasta_reads;
    auto fasta_reader1 = bioparser::createReader<FASTARead, bioparser::FastaReader>(read_file_path);
    fasta_reader1->read_objects(fasta_reads, static_cast<uint64_t>(-1));
    long number_of_reads = fasta_reads.size();
    printf("Reading file - Done.\n");

    unordered_map<uint64_t, vector<hashMinPair3>> lookup_map; // hash minimizera -> minimizeri svih sekvenci poredani po indeksu uzlazno
    std::vector<std::vector<minim>> mins_in_order(number_of_reads); // id sekvence -> poredani minimizeri sekvence po indeksu

    // create thread pool
    std::shared_ptr<thread_pool::ThreadPool> thread_pool_data;
    if(max_concurrent_threads < 0){
        thread_pool_data = thread_pool::createThreadPool();
    }else{
        thread_pool_data = thread_pool::createThreadPool(max_concurrent_threads);
    }
    // create storage for return values of find_overlaps_by_LIS
    std::vector<std::future<uint32_t >> thread_futures_data;

    printf("Colecting data [-]");
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    for (int i=0; i<number_of_reads; i++){
            //concurrency?
        thread_futures_data.emplace_back(thread_pool_data->submit_task(
                process_sequence4_id, fasta_reads[i]->get_data(),
                fasta_reads[i]->get_data_length(),
                i,
                w,
                k,
                std::ref(mins_in_order)
        ));
    }
    int i = 0;
    unordered_map<uint64_t ,uint32_t > min_occurences;
    //nogo planvector<uint64_t > nogos;
    uint32_t  number_of_minimizers=0;
    for (auto& it: thread_futures_data) {
        report_status("Collecting minimizers",i++, number_of_reads);
        it.wait();
        uint32_t id = it.get();
        number_of_minimizers+=add_to_lookup_table(id, mins_in_order[id], lookup_map, min_occurences);
    }
    //nogo plan reduce_minimizers(number_of_minimizers,min_occurences,nogos);
    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    printf("\rCollecting minimizers - Finished in %ld seconds.\n",chrono::duration_cast<chrono::seconds>( t2 - t1 ).count());

    sort_by_indices(lookup_map);
    chrono::high_resolution_clock::time_point t4 = chrono::high_resolution_clock::now();
    printf("Sorting finished in %ld seconds.\n",chrono::duration_cast<chrono::seconds>( t4 - t2 ).count());
    fflush(stdout);
    //nogo plan sort(nogos.begin(),nog  os.end());
    printf("Comparing sequences [-]");

    // create thread pool
    std::shared_ptr<thread_pool::ThreadPool> thread_pool_lis;
    if(max_concurrent_threads < 0){
        thread_pool_lis = thread_pool::createThreadPool();
    }else{
        thread_pool_lis = thread_pool::createThreadPool(max_concurrent_threads);
    }
    // create storage for return values of find_overlaps_by_LIS
    std::vector<std::future<void>> thread_futures_lis;


    FILE* output = fopen(result_file_path.c_str(),"w");
    chrono::high_resolution_clock::time_point t5 = chrono::high_resolution_clock::now();
    for (int i = 0; i < number_of_reads; ++i) {
        thread_futures_lis.emplace_back(thread_pool_lis->submit_task(
                lis_overlap_parallelization,
                i,
                std::ref(mins_in_order[i]),
                std::ref(lookup_map),
                4,
                std::ref(fasta_reads),
                output,
                std::ref(min_occurences)));
    }
    int j = 0;
    for (auto& it: thread_futures_lis) {
        report_status("Comparing sequences",j++, number_of_reads);
        it.wait();
    }

    chrono::high_resolution_clock::time_point t6 = chrono::high_resolution_clock::now();
    printf("\rComparing sequences - Finished in %ld seconds.\n", chrono::duration_cast<chrono::seconds>( t6 - t5 ).count());
    printf("Total execution time: %ld seconds\n", chrono::duration_cast<chrono::seconds>( t6 - t1 ).count());
    return 0;
}