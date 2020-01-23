#include <benchmark/benchmark.h>

#include "SearchStorageFixture.hpp"
#include <algorithm>
#include <iostream>

namespace {

template <class T>
constexpr void ignore(const T&) {}

using STLSearch = SearchStorageFixture<std::vector<A>>;

BENCHMARK_DEFINE_F(STLSearch, STLEqualRange)
(benchmark::State& state)
{
    size_t itr{}, sum{};
    for (auto _ : state) {
        ignore(_);
        auto id = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto first = std::equal_range(std::cbegin(records_), std::cend(records_), id, Comp{});
        if (first.first != first.second) {
            if ((*first.first).id_ != id) {
                std::cout << "error search" << std::endl;
            }
            benchmark::DoNotOptimize(sum += (*first.first).id_);
        } else {
            std::cout << "not found error" << std::endl;
        }
    }
}
/*
BENCHMARK_REGISTER_F(STLSearch, STLEqualRange)
        ->RangeMultiplier(0xF + 1)
        ->Ranges({{0xF + 1, 0xFFFFFF + 1}, {0, 2}})
        ->Complexity()
        ->MeasureProcessCPUTime()
        ->Threads(1);
        */
} // namespace