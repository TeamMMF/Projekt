//
// Created by filip on 06.11.17..
//

#include "gtest/gtest.h"
#include "Common.hpp"

namespace {

    TEST(ComplementTest, TestGuanine) {
        EXPECT_EQ('G', complement('C'));

    }

    TEST(ComplementTest, TestCytosine) {
        EXPECT_EQ('C', complement('G'));
    }

    TEST(ComplementTest, TestAdenine) {
        EXPECT_EQ('A', complement('T'));

    }

    TEST(ComplementTest, TestThymine) {
        EXPECT_EQ('T', complement('A'));

    }

    TEST(ComplementTest, TestInvalid) {
        EXPECT_ANY_THROW({complement('D');});

    }

    TEST(HitSortTest, TestSort) {
        cout << "tu";
       vector<minimizer_hit> vec;
        vec.push_back(make_tuple(0,0,1,0));
        vec.push_back(make_tuple(0,0,1,3));
        vec.push_back(make_tuple(2,1,2,5));
        vec.push_back(make_tuple(0,0,1,1));
        vec.push_back(make_tuple(0,0,1,4));
        vec.push_back(make_tuple(3,0,5,2));
        vec.push_back(make_tuple(0,3,1,0));
        vec.push_back(make_tuple(7,1,5,4));
        vec.push_back(make_tuple(3,1,5,4));
        vec.push_back(make_tuple(0,0,1,2));
        vec.push_back(make_tuple(0,0,0,1));
        vec.push_back(make_tuple(6,0,4,5));
        vec.push_back(make_tuple(0,1,0,10));
        vec.push_back(make_tuple(3,0,1,5));
        sort(vec.begin(),vec.end(),hit_comparator);
        for(auto t : vec){
            printf("(%d,%d,%d,%d)\n",
                   (int) get<0>(t),
                   (int) get<1>(t),
                   (int) get<2>(t),
                   (int) get<3>(t)
            );
        }
    }
}