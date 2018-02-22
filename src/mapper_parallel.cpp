//
// Created by filip on 19.12.17..
//

#include "ProjectConfig.h"
#include "Common.hpp"
#include "FASTARead.h"
#include "Dynamic.h"
#include <chrono>
#include "bioparser/bioparser.hpp"
#include <algorithm>
#include <include/thread_pool/thread_pool.hpp>

#define WINDOW_DEFAULT 5
#define KMER_DEFAULT 15
#define LIS_THRESHOLD 4

/**
 * A constant used for denoting progress to the user, simulates a rotating progress bar.
 */
const char *progress = "-\\|/";

/**
 * Shows the program usage instructions to the user.
 * @param arg The program's name.
 */
void show_usage(const string &arg) {
    fprintf(stdout, "%s Version %d.%d\n",
            arg.c_str(),
            SequenceOverlaping_VERSION_MAJOR,
            SequenceOverlaping_VERSION_MINOR
    );
    fprintf(stdout, "Usage: %s <option(s)> <readfile> <resultsfile>\n", arg.c_str());
    fprintf(stdout, "Available options are:\n");
    fprintf(stdout, "-h, --help\t\t\tshow usage instructions.\n");

}

/**
 * Adds all minimizers of one given sequence to a shared lookup map and
 * tracks how many times a given minimizer appears by using a shared counting map.
 * @param seq_id The id of the sequence to process.
 * @param minimizers A vector of all minimizers appearing in the sequence.
 * @param lookup_map A lookup map shared by all sequences, provides a vector of all minimizers for a particular hash.
 * @param occurences A shared map used for counting minimizer occurrences.
 */
void add_to_lookup_map(uint32_t seq_id,
                       std::vector<minim> &minimizers,
                       std::unordered_map<uint64_t, vector<hashMinPair3>> &lookup_map,
                       std::unordered_map<uint64_t, uint32_t> &occurences) {
    uint64_t len = minimizers.size();
    for (uint32_t j = 0; j < len; j++) {
        minim min = minimizers[j];
        auto it = lookup_map.find(min.hash);

        if (it == lookup_map.end()) {
            std::vector<hashMinPair3> vec;
            vec.emplace_back((hashMinPair3) {seq_id, min.index});
            lookup_map.emplace(min.hash, vec);
        } else {
            it->second.emplace_back((hashMinPair3) {seq_id, min.index});
        }

        occurences[min.hash]++;
    }
}

/**
 * Reports a status of a given operation to the user using a rotating progress bar.
 * @param operation The operation being executed.
 * @param curr The current number of processed elements.
 * @param total The total number of elements to be processed.
 */
void report_status(const char *operation, int curr, long total) {
    long ratio = 100 * curr / total;
    fprintf(stdout, "\r%s [%c] %ld%c", operation, progress[curr % 4], ratio, '%');
    fflush(stdout);
}

/**
 * Finds all overlaps for a sequence and outputs them to a specified output stream.
 * @param query_id The query sequence's id.
 * @param minimizer_hashes A vector of all minimizers appearing in the query sequence.
 * @param lookup_map A map mapping a hash to all corresponding minimizers (from all sequences)
 * @param fasta_reads A vector of FASTARead objects containing data required by the PAF format
 * @param output The result output stream. This is where all of the found matches are printed (in PAF format)
 * @param occurences A map denoting the number of occurences for a given minimizer hash (based on all sequences)
 */
void find_overlaps(int query_id,
                   vector<minim> &minimizer_hashes,
                   unordered_map<uint64_t, vector<hashMinPair3>> &lookup_map,
                   vector<unique_ptr<FASTARead>> &fasta_reads,
                   FILE *output,
                   unordered_map<uint64_t, uint32_t> &occurences) {

    vector<pair<int, bool>> result = find_overlaps_by_LIS(query_id, minimizer_hashes, lookup_map,
                                                          LIS_THRESHOLD, occurences);
    for (auto res : result) {
        fprintf(output, "%s\t%d\t%d\t%d\t%c\t%s\t%d\n",
                fasta_reads[query_id]->get_name(),
                fasta_reads[query_id]->get_data_length(),
                0,
                0,
                res.second ? '+' : '-',
                fasta_reads[res.first]->get_name(),
                fasta_reads[res.first]->get_data_length()
        );
    }
}

/**
 * Starts the program.
 * @param argc The number of command line arguments.
 * @param argv The command line arguments
 * @return 0 if the execution succeeds, 1 otherwise
 */
int main(int argc, char const *argv[]) {

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
    if (max_concurrent_threads < 0) {
        printf("Running on an unlimited number of concurrent threads.\n");
    } else {
        printf("Running on %d concurrent threads.\n", max_concurrent_threads);
    }

    // Reads the file and creates a FASTARead object for every read
    vector<unique_ptr<FASTARead>> fasta_reads;
    auto fasta_reader1 = bioparser::createReader<FASTARead, bioparser::FastaReader>(read_file_path);
    fasta_reader1->read_objects(fasta_reads, static_cast<uint64_t>(-1));
    unsigned long number_of_reads = fasta_reads.size();

    printf("Reading file - Done.\n");

    // maps a hash to a vector of corresponding minimizers (taken from all available sequences)
    unordered_map<uint64_t, vector<hashMinPair3>> lookup_map;
    // maps a sequence id to an ascenging list of its minimizers
    std::vector<std::vector<minim>> ordered_minimizers(
            number_of_reads); // id sekvence -> poredani minimizeri sekvence po indeksu

    // create thread pool for collecting minimizers
    std::shared_ptr<thread_pool::ThreadPool> thread_pool_data;
    if (max_concurrent_threads < 0) {
        thread_pool_data = thread_pool::createThreadPool();
    } else {
        thread_pool_data = thread_pool::createThreadPool(max_concurrent_threads);
    }
    std::vector<std::future<uint32_t >> thread_futures_data;

    printf("Colecting data [-]");
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    for (int i = 0; i < number_of_reads; i++) {
        thread_futures_data.emplace_back(thread_pool_data->submit_task(
                process_sequence4_id, fasta_reads[i]->get_data(),
                fasta_reads[i]->get_data_length(),
                i,
                WINDOW_DEFAULT,
                KMER_DEFAULT,
                std::ref(ordered_minimizers)
        ));
    }

    unordered_map<uint64_t, uint32_t> minimizer_occurences;
    int coll_progress = 0;
    for (auto &it: thread_futures_data) {
        report_status("Collecting minimizers", coll_progress++, number_of_reads);
        it.wait();
        uint32_t id = it.get();
        add_to_lookup_map(id, ordered_minimizers[id], lookup_map, minimizer_occurences);
    }

    chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
    printf("\rCollecting minimizers - Finished in %ld seconds.\n",
           chrono::duration_cast<chrono::seconds>(t2 - t1).count());

    sort_by_indices(lookup_map);
    chrono::high_resolution_clock::time_point t4 = chrono::high_resolution_clock::now();
    printf("Sorting finished in %ld seconds.\n", chrono::duration_cast<chrono::seconds>(t4 - t2).count());
    fflush(stdout);

    printf("Comparing sequences [-]");

    // create thread pool for comparing sequences
    std::shared_ptr<thread_pool::ThreadPool> thread_pool_lis;
    if (max_concurrent_threads < 0) {
        thread_pool_lis = thread_pool::createThreadPool();
    } else {
        thread_pool_lis = thread_pool::createThreadPool(max_concurrent_threads);
    }
    std::vector<std::future<void>> thread_futures_lis;


    FILE *output = fopen(result_file_path.c_str(), "w");
    chrono::high_resolution_clock::time_point t5 = chrono::high_resolution_clock::now();
    for (int i = 0; i < number_of_reads; ++i) {
        thread_futures_lis.emplace_back(thread_pool_lis->submit_task(
                find_overlaps,
                i,
                std::ref(ordered_minimizers[i]),
                std::ref(lookup_map),
                std::ref(fasta_reads),
                output,
                std::ref(minimizer_occurences)));
    }

    int cmp_progress = 0;
    for (auto &it: thread_futures_lis) {
        report_status("Comparing sequences", cmp_progress++, number_of_reads);
        it.wait();
    }

    chrono::high_resolution_clock::time_point t6 = chrono::high_resolution_clock::now();
    printf("\rComparing sequences - Finished in %ld seconds.\n",
           chrono::duration_cast<chrono::seconds>(t6 - t5).count());
    printf("Total execution time: %ld seconds\n", chrono::duration_cast<chrono::seconds>(t6 - t1).count());

    return 0;
}