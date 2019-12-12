#include "gtest/gtest.h"

namespace my {

namespace testing {

    namespace {

        TEST(FactorialTest, Negative)
        {

            EXPECT_EQ(1, 1);
        }

        TEST(FactorialTest, Positive)
        {

            EXPECT_EQ(4, 2);
        }

        TEST(TestMainTest, ShouldSucceed)
        {
        }

    }
}
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
