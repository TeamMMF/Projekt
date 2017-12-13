//
// Created by filip on 12.12.17..
//

#ifndef SEQUENCEOVERLAPING_FASTAREAD_H
#define SEQUENCEOVERLAPING_FASTAREAD_H

#include <cstdlib>
#include <cstring>
#include <cstdint>

class FASTARead {

    char *name;
    char *data;

    uint32_t name_length;
    uint32_t data_length;

public:
    FASTARead(uint64_t object_id,
              const char *name,
              uint32_t name_length,
              const char *data,
              uint32_t data_length) {

        this->name_length = name_length;

        this->name = (char *) malloc((name_length + 1) * sizeof(char));
        strncpy(this->name, name, name_length);
        this->name[name_length] = '\0';

        this->data_length = data_length;
        this->data = (char *) malloc((data_length + 1) * sizeof(char));
        strncpy(this->data, data, data_length);
        this->data[data_length] = '\0';
    }

    const char *get_data();

    const char *get_name();

    uint32_t get_data_length();


};


#endif //SEQUENCEOVERLAPING_FASTAREAD_H
