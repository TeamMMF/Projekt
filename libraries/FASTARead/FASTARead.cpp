//
// Created by filip on 12.12.17..
//

#include "FASTARead.h"

const char* FASTARead::get_name(){
    return name_.c_str();
};

const char* FASTARead::get_data(){
    return data_.c_str();
};

uint32_t FASTARead::get_data_length(){
    return data_.size();
}

uint32_t FASTARead::get_name_length() {
    return name_.size();
}
