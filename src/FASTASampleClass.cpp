//
// Created by filip on 06.11.17..
//

#include <string.h>
#include <string>
#include <sstream>

// define a class for sequences in FASTA format
class FASTASampleClass {

    char *name;
    char *data;

    uint32_t name_length;
    uint32_t data_length;

public:
    FASTASampleClass(uint64_t object_id,
                     const char *name,
                     uint32_t name_length,
                     const char *data,
                     uint32_t data_length) {

        this->name_length = name_length;
        this->name = new char[name_length];
        strncpy(this->name, name, name_length);

        this->data_length = data_length;
        this->data = new char[data_length];
        strncpy(this->data, data, data_length);
    }

    std::string get_description();
    std::string get_data();
};

std::string FASTASampleClass::get_description() {
    std::ostringstream oss;

    oss << "Name: ";
    for (int i = 0; i < name_length; ++i) {
        oss << name[i];
    }

    oss << '\n';
    oss << "Data: ";
    for (int i = 0; i < data_length; ++i) {
        oss << data[i];
    }

    return oss.str();
}

std::string FASTASampleClass::get_data() {
    std::ostringstream oss;

    for (int i = 0; i < data_length; ++i) {
        oss << data[i];
    }

    return oss.str();
};