//
// Created by filip on 12.12.17..
//

#include "FASTARead.h"

const char* FASTARead::get_name(){
    return name;
};

const char* FASTARead::get_data(){
    return data;
};

uint32_t FASTARead::get_data_length(){
    return data_length;
}

uint32_t FASTARead::get_name_length() {
    return name_length;
}