// Performance.cpp : Defines the entry point for the application.
//

#include "Performance.h"
#include <benchmark/benchmark.h>
#include <random>
#include <vector>

struct A {
    A(uint32_t id)
        : id_(id)
    {
    }
    uint32_t id_;
};

size_t bsearch(const std::vector<A*>& resource, uint32_t key)
{
    size_t low = 0;
    size_t high = std::min<size_t>(key, resource.size());
    size_t mid = high / 2;

    while (low < high) {
        auto val = resource[mid]->id_;

        if (key < val) {
            high = mid;
        } else {
            if (key == val) {
                break;
            }
            low = mid + 1;
        }

        mid = low + (high - low) / 2;
    }
    return resource.size() && resource[mid]->id_ == key ? mid : resource.size();
}

size_t bsearch_bf(const std::vector<A*>& resource, uint32_t key)
{
    size_t lower = 0;
    size_t n = std::min<size_t>(resource.size(), key);

    while (auto half = n / 2) {
        auto middle = lower + half;
        auto val = resource[middle]->id_;
        lower = (val <= key) ? middle : lower;
        n -= half;
    }
    return (resource[lower]->id_ == key) ? lower : resource.size();
}

static void BM_SomeFunction(benchmark::State& state)
{
    size_t test_cases = 0xFFFF;
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist6(1, test_cases); // distribution in range [1, 6]
    std::vector<A*> vec(test_cases, nullptr);
    std::cout << vec.size() << std::endl;
    for (size_t i = 0; i < vec.size(); i++) {
        vec[i] = new A(static_cast<uint32_t>(i));
    }
    // Perform setup here
    for (auto _ : state) {
        // This code gets timed
        //auto i = std::lower_bound(vec.cbegin(), vec.cend(), dist6(rng)); 
        auto i = bsearch(vec, dist6(rng));
        //
    }
}
// Register the function as a benchmark
BENCHMARK(BM_SomeFunction);
// Run the benchmark
BENCHMARK_MAIN();