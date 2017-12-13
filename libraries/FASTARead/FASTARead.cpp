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

int FASTARead::get_data_length(){
    return data_length;
}
