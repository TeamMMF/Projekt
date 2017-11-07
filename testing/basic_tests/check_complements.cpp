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
}