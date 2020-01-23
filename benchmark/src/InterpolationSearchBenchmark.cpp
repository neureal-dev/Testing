#include <benchmark/benchmark.h>

#include <array>
#include <iostream>
#include <memory>
#include <vector>

#include "InterpolationSearch.hpp"
#include "SearchStorageFixture.hpp"

namespace {

template <class T>
constexpr void ignore(const T&) {}

using InterpolationSearch = SearchStorageFixture<std::vector<A>>;

BENCHMARK_DEFINE_F(InterpolationSearch, Hybrid)
(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        ignore(_);
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto low = mtl::HybridInterpolationSearch(
                std::cbegin(records_),
                std::cend(records_),
                rid,
                Comp{},
                //            [](const A* first, const A* last, uint32_t key) -> double
                //            { return (double(key) - first->GetId()) / (last->GetId() -
                //            first->GetId()); });
                [](const A& first, const A& last, uint32_t key) -> double {
                    return (double(key) - first.id_) / (last.id_ - first.id_);
                });

        if (low != std::cend(records_)) {
            if ((*low).id_ != rid) {
                std::cout << "error search" << std::endl;
            }
            benchmark::DoNotOptimize(sum += (*low).id_);
        } else {
            std::cout << "not found error" << std::endl;
        }
    }
}
/*BENCHMARK_REGISTER_F(InterpolationSearch, Hybrid)
        ->RangeMultiplier(0xF + 1)
        ->Ranges({{0xF + 1, 0xFFFFFF + 1}, {0, 2}})
        ->Complexity()
        ->MeasureProcessCPUTime()
        ->Threads(6);
*/
BENCHMARK_DEFINE_F(InterpolationSearch, Pure)
(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        ignore(_);
        /*
        int s = 0;
        for (auto itr = records_.begin(); itr != records_.end(); ++itr) {
            benchmark::DoNotOptimize(s += itr->id_);
        }
        //std::cout << s << std::endl;
        */
        
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto low = mtl::InterpolationSearch(
                std::cbegin(records_),
                std::cend(records_),
                rid,
                Comp{},
                //            [](const A* first, const A* last, uint32_t key) -> double
                //            { return (double(key) - first->GetId()) / (last->GetId() -
                //            first->GetId()); });
                [](const A& first, const A& last, uint32_t key) -> double {
                    return (double(key) - first.id_) / (last.id_ - first.id_);
                });

        if (low != std::cend(records_)) {
            if ((*low).id_ != rid) {
                std::cout << "error search" << std::endl;
            }
            benchmark::DoNotOptimize(sum += (*low).id_);
        } else {
            std::cout << "not found error" << std::endl;
        }
        
    }
}

BENCHMARK_REGISTER_F(InterpolationSearch, Pure)
        ->RangeMultiplier(0xF + 1)
        ->Ranges({{0xFF + 1, 0xFFFFFF + 1}, {0, 2}})
        ->Complexity()
        ->MeasureProcessCPUTime()
        ->Threads(1);
        
} // namespace
