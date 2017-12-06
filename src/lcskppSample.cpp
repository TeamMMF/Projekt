//
// Created by filip on 06.12.17..
//
#include <iostream>
#include "lcskpp.h"

using namespace std;

int main(){
    int length;
    lcskpp("ABXXXCDE","ABYYYCDE",2,&length,NULL);
    cout <<length<<endl;
    return 0;
}

