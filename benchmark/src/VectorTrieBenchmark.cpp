#include <benchmark/benchmark.h>

#include <iostream>

#include "SearchStorageFixture.hpp"
#include "VectorTrie.hpp"

namespace {

template <class T>
constexpr void ignore(const T&) {}

using SearchOptimization = SearchStorageFixture<mtl::PersistentVector<A, 17, 0>>;

BENCHMARK_DEFINE_F(SearchOptimization, RBTrie)
(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        ignore(_);
        int s = 0;
        for (auto itr = records_.leafs_.begin(); itr != records_.leafs_.end(); ++itr) {
           // std::cout << (*itr).id_ << std::endl;
            if (*itr)
            s += (*itr)->id_;
            //s +=;
        }
        /*
        for (auto i = 0; i < std::numeric_limits<int>::max(); i++) {
            auto j = records_.getNode(i);
            if (j)
            s += j->id_;
        }
        */
        //std::cout << s << std::endl;
        /*

        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto rec = records_.getNode(rid);
        if (rec) {
            benchmark::DoNotOptimize(sum += rec->id_);
        } else {
            std::cout << "record not found error" << std::endl;
        }
        */
    }
}

BENCHMARK_REGISTER_F(SearchOptimization, RBTrie)
        ->RangeMultiplier(0xFFFF + 1)
        ->Ranges({{0xFF + 1, 0xFF + 1}, {0, 2}})
        ->Complexity()
        ->MeasureProcessCPUTime()
        ->Threads(1);

} // namespace
