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

int HybridSearch(const std::vector<A*>& a, uint32_t x)
{
    int64_t left = 0;
    int64_t right = a.size() - 2;
    int64_t Inter, Mid;

    while (left < right) { //1 Comps
        Inter = left + (((static_cast<int64_t>(x) - a[left]->id_) / (a[right]->id_ - a[left]->id_)) * (right - left));
        if (x > a[Inter]->id_) { //2Comps
            Mid = (Inter + right) / 2;
            if (x <= a[Mid]->id_) //3Comps
            {
                left = Inter + 1;
                right = Mid;
            } else {
                left = Mid + 1;
            }
        } else if (x < a[Inter]->id_) { //3Comps
            Mid = (Inter + left) / 2;
            if (x >= a[Mid]->id_) { //4Comps
                left = Mid;
                right = Inter - 1;
            } else {
                right = Mid - 1;
            }
        } else {
            return Inter;
        }
    } //endwhile
    if (left < a.size() && x == a[left]->id_)
        return left;
    return -1; //notfound
}

size_t bsearch(const std::vector<A*>& resource, uint32_t key)
{
    size_t low = 0;
    size_t high = resource.size();
    if (high > 1) {
        high = high - 1;
        //high = std ::min<size_t>(key, high - 1)
        //high = std ::min<size_t>(1ull + key - resource[low]->id_, high - 1);
    }
    //size_t high = resource.size() - 1;
    size_t ind = -1;
    size_t mid = high >> 1;

    while (low < high) {
        //auto mid = low + ((high - low) >> 1);
        auto val = resource[mid]->id_;

        if (key < val) {
            high = mid;
        } else {
            if (key == val) {
                ind = mid;
                break;
            }
            low = mid + 1;
        }

        mid = low + ((high - low) >> 1);
    }
    return ind;
    //return resource.size() && resource[mid]->id_ == key ? mid : resource.size();
}

size_t bsearch_bf(const std::vector<A*>& resource, uint32_t key)
{
    size_t lower = 0;
    size_t n = resource.size();
    uint32_t val = 0;
    while (auto half = (n >> 1)) {
        auto middle = lower + half;
        val = resource[middle]->id_;
        lower = (val > key) ? lower : middle;
        n -= half;
    }
    return (val == key) ? lower : resource.size();
}

static void BM_SomeFunction(benchmark::State& state)
{
    size_t test_cases = 0xFFF;
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist6(1ull, test_cases); // distribution in range [1, 6]
    std::vector<A*> vec(test_cases, nullptr);
    std::cout << vec.size() << std::endl;
    for (size_t i = 0; i < vec.size(); i++) {
        vec[i] = new A(static_cast<uint32_t>(i * 2));
    }
    // Perform setup here
    for (auto _ : state) {
        // This code gets timed
        //auto i = std::lower_bound(vec.cbegin(), vec.cend(), dist6(rng)); 
        auto i = bsearch_bf(vec, dist6(rng) * 2);
        //auto i = HybridSearch(vec, dist6(rng) + 65536);
        //
    }
}
// Register the function as a benchmark
BENCHMARK(BM_SomeFunction);
// Run the benchmark
BENCHMARK_MAIN();