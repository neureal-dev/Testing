#include <benchmark/benchmark.h>

#include <iostream>
#include <vector>
#include <array>
#include <memory>

#include "InterpolationSearch.hpp"
#include "SearchStorageFixture.hpp"

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

struct Comp {

    inline bool operator()(const A* s, uint32_t i) const noexcept { return s->GetId() < i; }

    inline bool operator()(uint32_t i, const A* s) const noexcept { return i < s->GetId(); }

    inline bool operator()(const A& s, uint32_t i) const noexcept { return s.GetId() < i; }

    inline bool operator()(uint32_t i, const A& s) const noexcept { return i < s.GetId(); }

    inline bool operator()(const std::unique_ptr<A>& s, uint32_t i) const noexcept { return s->GetId() < i; }

    inline bool operator()(uint32_t i, const std::unique_ptr<A>& s) const noexcept { return i < s->GetId(); }
};

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
        auto low = mtl::HybridInterpolationSearch(std::cbegin(records_),
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
BENCHMARK_REGISTER_F(InterpolationSearch, Hybrid)
        ->RangeMultiplier(0xF + 1)
        ->Ranges({{0xF + 1, 0xFFFFFF + 1}, {0, 2}})
        ->Complexity()
        //->MeasureProcessCPUTime()
        ->Threads(6);

} // namespace
