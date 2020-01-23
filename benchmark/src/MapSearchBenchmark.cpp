#include <benchmark/benchmark.h>

#include <iostream>
#include <map>
#include <unordered_map>

#include "SearchStorageFixture.hpp"
#include "VectorTrie.hpp"

namespace {

template <typename T>
struct bmmap {
    using value_type = T;
    using pointer_type = typename std::add_pointer<T>::type;

    void push_back(const value_type& value)
    {
        map.emplace(value.id_, value);
    }

    pointer_type getNode(uint32_t id)
    {
        return &(map.find(id)->second);
    }

    void clear()
    {
        map.clear();
    }

    void shrink_to_fit()
    {
        //
    }

    std::unordered_map<uint32_t, value_type> map;
};

template <class T>
constexpr void ignore(const T&) {}

using SearchOptimization = SearchStorageFixture<bmmap<A>>;

BENCHMARK_DEFINE_F(SearchOptimization, STLMap)
(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        ignore(_);
        auto rid = searches_[(searches_.size() - ++itr) % searches_.size()];
        auto rec = records_.getNode(rid);
        if (rec) {
            benchmark::DoNotOptimize(sum += rec->id_);
        } else {
            std::cout << "record not found error" << std::endl;
        }
    }
}
/*
BENCHMARK_REGISTER_F(SearchOptimization, STLMap)
        ->RangeMultiplier(0xF + 1)
        ->Ranges({{0xF + 1, 0xFFFFFF + 1}, {0, 2}})
        ->Complexity()
        ->MeasureProcessCPUTime()
        ->Threads(1);
*/
} // namespace
