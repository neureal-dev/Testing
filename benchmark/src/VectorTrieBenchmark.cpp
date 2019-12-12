#include <benchmark/benchmark.h>

#include <iostream>

#include "SearchStorageFixture.hpp"
#include "VectorTrie.hpp"

namespace {

struct A {
    explicit A(uint32_t id)
            : id_(id), id1_{}, id2_{}, id4_{}
    {
    }

    [[nodiscard]] uint32_t GetId() const noexcept { return id_; }

    uint32_t id_;
    uint64_t id1_;
    uint64_t id2_;
    std::array<uint64_t, 3> id4_;
};

template <class T>
constexpr void ignore(const T&) {}

using SearchOptimization = SearchStorageFixture<mtl::PersistentVector<A, 2, 16>>;

BENCHMARK_DEFINE_F(SearchOptimization, RBTrie)
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

BENCHMARK_REGISTER_F(SearchOptimization, RBTrie)
        ->RangeMultiplier(0xF + 1)
        ->Ranges({{0xF + 1, 0xFFFFFF + 1}, {0, 2}})
        ->Complexity()
        //->MeasureProcessCPUTime()
        ->Threads(6);

} // namespace
