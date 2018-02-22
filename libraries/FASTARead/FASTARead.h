//
// Created by filip on 12.12.17..
//

#ifndef SEQUENCEOVERLAPING_FASTAREAD_H
#define SEQUENCEOVERLAPING_FASTAREAD_H

#include <cstdlib>
#include <string>
#include <cstdint>

/**
 * The class FASTARead encapsulates information about a read taken from a FASTA file.
 * The class with this format was required by the library bioparser (https://github.com/rvaser/bioparser) used for reading data from FASTA files.
 */
class FASTARead {
    uint64_t id_;
    std::string name_;
    std::string data_;

public:
    FASTARead(uint64_t id,
              const char *name,
              uint32_t name_length,
              const char *data,
              uint32_t data_length) : id_(id), name_(name, name_length), data_(data, data_length) {
    }

    const char*  get_data();

    const char*  get_name();

    uint32_t get_data_length();

    uint32_t get_name_length();

};


#endif //SEQUENCEOVERLAPING_FASTAREAD_H
