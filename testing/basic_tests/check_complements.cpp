//
// Created by filip on 06.11.17..
//

#include "gtest/gtest.h"
#include "Common.h"

namespace {

    TEST(CommonTest, TestGuanine) {
        EXPECT_EQ('G', complement('C'));

    }

    TEST(CommonTest, TestCitosine) {
        EXPECT_EQ('C', complement('G'));
    }

    TEST(CommonTest, TestAdenin) {
        EXPECT_EQ('A', complement('T'));

    }

    TEST(CommonTest, TestTimine) {
        EXPECT_EQ('T', complement('A'));

    }
}