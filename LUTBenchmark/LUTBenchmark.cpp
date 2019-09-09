// Test.cpp : Defines the entry point for the application.
//

#include "LUTBenchmark.h"
#include <benchmark/benchmark.h>
#include <algorithm>
#include <chrono>
#include <mutex>
#include <random>
#include <shared_mutex>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

class A {
public:
    A() : id_(0) {}
    A(uint32_t id) : id_(id) {}
    uint32_t getId() const { return id_; }

    uint32_t id_;
    uint64_t id1_;
    uint64_t id2_;
    uint64_t id3_;
    uint64_t id4_;
};

class VectorSearchFixture : public benchmark::Fixture {
public:
    void SetUp(const ::benchmark::State& state)
    {
        records_ = std::vector<A>(state.range_x());
        searches_ = std::vector<uint32_t>(state.range_x(), 0);
        benchmark::DoNotOptimize(records_.data());
        benchmark::DoNotOptimize(searches_.data());
        std::mt19937 rng;
        rng.seed(std::random_device()());
        std::uniform_int_distribution<uint32_t> uniform_dist(0, state.range_x());

        for (int64_t i = state.range_x() - 1; i >= 0; --i) {
            records_[i] = (i + 65536);
            searches_[i] = uniform_dist(rng);
        }
    }

    void TearDown(const ::benchmark::State& state)
    {
        for (auto item : records_) {
            //delete item;
        }
        records_.clear();
    }

//    private:
    std::vector<A> records_;
    std::vector<uint32_t> searches_;
};

enum locate_t { EQUAL, LEFT, RIGHT };

int nary_search(const std::vector<A>& a, int const key, int const N)
{
    std::vector<int> mid(N + 1);
    std::vector<locate_t> locate(N + 2);

    locate[0] = RIGHT;
    locate[N + 1] = LEFT;

    int lo = 0;
    int hi = a.size() - 1;
    int pos = -1;

    while (lo <= hi && pos == -1) {
        mid[0] = lo - 1;

        double const step = (hi - lo + 1) / (N + 1);

        for (int i = 1; i <= N; i++) {
            int const offset = step * i + (i - 1);
            int const lmid = mid[i] = lo + static_cast<int>(offset);

            if (lmid <= hi) {
                if (a[lmid].id_ > key) {
                    locate[i] = LEFT;
                } else if (a[lmid].id_ < key) {
                    locate[i] = RIGHT;
                } else {
                    locate[i] = EQUAL;
                    pos = lmid;
                }
            } else {
                mid[i] = hi + 1;
                locate[i] = LEFT;
            }
        }
        for (int i = 1; i <= N; i++) {
            if (locate[i] != locate[i - 1]) {
                lo = mid[i - 1] + 1;
                hi = mid[i] - 1;
            }
        }
        if (locate[N] != locate[N + 1]) {
            lo = mid[N] + 1;
        }
    }
    return pos;
}


int binary_cmov(const std::vector<A>& records, uint32_t rid)
{
    size_t min = 0, max = records.size();
    while (min < max) {
        size_t middle = (min + max) >> 1;
        //middle /= 2;
        auto val = records[middle].id_;
        min = rid > val ? middle + 1 : min;
        max = rid > val ? max : middle;
    }
    return min;
}

int binary_cmov2(std::vector<A>& const records, uint32_t const rid)
{
    size_t min = 0, max = records.size();
    while (min < max) {
        size_t middle = (min + max) >> 1;
        size_t middle1 = middle + 1;
        min = rid > records[middle].id_ ? middle1 : min;
        max = rid <= records[middle].id_ ? middle : max;
    }
    return min;
}

size_t b2search(std::vector<A>& const records, uint32_t const rid)
{
    size_t size = records.size();
    size_t low = 0;
    //size_t middle;
    while (size_t half = size >> 1) {
        size_t middle = low + half;
        auto val = records[middle].id_;
        low = (val <= rid) ? middle : low;
        size = size - half;
    }
    return (size && records[low].id_ == rid) ? low : -1;
}

inline size_t bb2search(std::vector<A> const& records, uint32_t const rid, size_t low, size_t high)
{
    while (size_t half = high >> 1) {
        size_t middle = low + half;
        auto val = records[middle].id_;
        low = (val <= rid) ? middle : low;
        high = high - half;
    }
    return (high && records[low].id_ == rid) ? low : -1;
}


int lsearch(const std::vector<A>& records, uint32_t key)
{
    for (size_t i = 0; i < records.size(); i++) {
        if (records[i].id_ == key) {
            return i;
        }
    }
    return -1;
}

inline size_t bbsearch(std::vector<A> const& resource, uint32_t key, size_t low, size_t high)
{
    size_t ind = -1;
    size_t mid = high >> 1;

    while (low < high) {
        //size_t mid = (low + high) >> 1;
        auto val = resource[mid].id_;

        if (key > val) {
            low = mid + 1;
        } else {
            if (key == val) {
                ind = mid;
                break;
            }
            high = mid;
        }
        mid = (low + high) >> 1; 
    }
    return ind;
}

size_t bsearch(const std::vector<A>& data, uint32_t const key)
{
    size_t low = 0;
    
	size_t high = data.size();
    size_t ind = -1;
    size_t mid = high >> 1;

    while (low < high) {
        const auto val = data[mid].id_;
        if (key > val) {
            low = mid + 1;
        } else {
            if (key == val) {
                ind = mid;
                break;
            }
            high = mid;
        }
        mid = (low + high) >> 1;
    }
    return ind;
}


// size_t bsearch(const std::vector<A>& data, uint32_t const key) { return bbsearch(data, key, 0, data.size()); }

size_t fallback_bsearch(const std::vector<A>& data, uint32_t key)
{
    size_t low = 0;
    size_t high = data.size();
    size_t mid = high >> 1;
    // size_t ind = -1;
    size_t size = high;
    while (size > 0xFF) {
        // auto mid = (low + high) >> 1;
        auto val = data[mid].id_;
        if (key > val) {
            low = mid + 1;
        } else {
            if (key == val) {
                return mid;
                // low = high = mid;
                // break;
            }
            high = mid;
        }
        mid = (low + high) >> 1;
        size = high - low;
    }
    return bbsearch(data, key, low, high);
}

size_t combine_search(std::vector<A> const& data, uint32_t const rid)
{
    return (data.size() > 0xFFFFFull) ? bbsearch(data, rid, 0ull, data.size()) : bb2search(data, rid, 0ull, data.size());
}

size_t b2search_original(const std::vector<A>& records, uint32_t rid)
{
    int64_t numrecs = records.size();
    int64_t foundindex = -1;
    int64_t a = 0;
    int64_t z = numrecs - 1;
    int64_t k;
    uint32_t id;

    if (numrecs == 0)
        return records.size();
    while (a < z) {
        k = a + ((z - a) >> 1);

        id = records[k].id_;

        if (id > rid) {
            z = k;
        } else if (id < rid) {
            a = k + 1;
        } else {
            foundindex = k;
            break;
        }
        
           if (z - a <= 1) {
               if (a == 0 && records[a].id_ == rid)
                   foundindex = a;
               if (z == numrecs - 1 && records[z].id_ == rid)
                   foundindex = z;
               break;
           }
          
    }

    return foundindex != -1 ? foundindex : records.size();
}

/*
BENCHMARK_DEFINE_F(VectorSearchFixture, lookup_table)(benchmark::State& state)
{
    size_t itr{}, sum{};
    for (auto _ : state) {
        auto rid = searches_[itr++ & state.range_x() - 1];
        benchmark::DoNotOptimize(sum += records_[rid].id_);
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, lookup_table)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1);
*/

BENCHMARK_DEFINE_F(VectorSearchFixture, binary_cmov)(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        auto rid = searches_[itr++ & state.range_x() - 1];
        auto low = binary_cmov(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low].id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, binary_cmov)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1);


BENCHMARK_DEFINE_F(VectorSearchFixture, nary_search)(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        auto rid = searches_[itr++ & state.range_x() - 1];
        auto low = nary_search(records_, rid, 3);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low].id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, nary_search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1);


BENCHMARK_DEFINE_F(VectorSearchFixture, lsearch)(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        auto rid = searches_[itr++ & state.range_x() - 1];
        auto low = lsearch(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low].id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, lsearch)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1);

BENCHMARK_DEFINE_F(VectorSearchFixture, binary_fallback)(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        auto rid = searches_[itr++ & state.range_x() - 1];
        auto low = fallback_bsearch(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low].id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, binary_fallback)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1);

BENCHMARK_DEFINE_F(VectorSearchFixture, binary_cmov2)(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        auto rid = searches_[itr++ & state.range_x() - 1];
        auto low = binary_cmov2(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low].id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, binary_cmov2)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1);

BENCHMARK_DEFINE_F(VectorSearchFixture, bsearch)(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        auto rid = searches_[itr++ & state.range_x() - 1];
        auto low = bsearch(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low].id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, bsearch)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1);

BENCHMARK_DEFINE_F(VectorSearchFixture, combine_search)(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        auto rid = searches_[itr++ & state.range_x() - 1];
        auto low = combine_search(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low].id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, combine_search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1);

BENCHMARK_DEFINE_F(VectorSearchFixture, b2search)(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        auto rid = searches_[itr++ & state.range_x() - 1];
        auto low = b2search(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low].id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, b2search)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1);

BENCHMARK_DEFINE_F(VectorSearchFixture, b2search_original)(benchmark::State& state)
{
    uint64_t sum{}, itr{};
    for (auto _ : state) {
        auto rid = searches_[itr++ & state.range_x() - 1];
        auto low = b2search_original(records_, rid);
        if (low < records_.size()) {
            benchmark::DoNotOptimize(sum += records_[low].id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, b2search_original)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1);

/*
BENCHMARK_F(VectorSearchFixture, test_vector_bsearch)(benchmark::State& state)
{
    uint64_t sum{};
    size_t k, a = 0, z = records_.size(), foundindex, numrecs = records_.size(), itr = 0;
    for (auto _ : state) {
        auto rid = searches_[itr++ & state.range_x() - 1];
        // auto first = std::lower_bound(std::begin(vect), std::end(vect), id, [](const A lhs, auto i) . bool { return
        // lhs.id_ < i; });
        a = 0;
        z = records_.size();
        while (z - a > 1) {
            k = (a + z) >> 1;

            uint32_t id = records_[k].id_;
            if (id < rid) {
                a = k;
            } else if (id > rid) {
                z = k;
            } else {
                foundindex = k;
                break;
            }
        }

        if (z - a <= 1) {
            if (a == 0 && records_[a].id_ == rid)
                foundindex = a;
            if (z == numrecs - 1 && records_[z].id_ == rid)
                foundindex = z;
            //    break;
        }

        sum += records_[foundindex].id_;
        //  }
        // std::cout << sum << std::endl;
    }
}
*/
BENCHMARK_DEFINE_F(VectorSearchFixture, std_lower_bound)(benchmark::State& state)
{
    size_t itr{}, sum{};
    for (auto _ : state) {
        auto id = searches_[itr++ & state.range_x() - 1];
        auto first = std::lower_bound(records_.begin(), records_.end(), id,
                                      [](const A& lhs, auto i) -> bool { return lhs.id_ < i; });
        if (first != records_.end()) {
            benchmark::DoNotOptimize(sum += (*first).id_);
        }
    }
}
BENCHMARK_REGISTER_F(VectorSearchFixture, std_lower_bound)->RangeMultiplier(0xF + 1)->Range(0xF + 1, 0xFFFFFF + 1);
/*
void test_umap()
{
    const size_t test_size = 1000000;
    std::unordered_map<uint32_t, A> umap;
    umap.reserve(test_size);
    for (size_t i = 0; i < test_size; i++) {
        umap[static_cast<uint32_t>(i)] = new A(static_cast<uint32_t>(i));
    }
    std::random_device r;
    std::default_random_engine e1(r());
    std::uniform_int_distribution<int> uniform_dist(1, test_size - 1);
    uint64_t sum = 0;
    for (size_t i = 0; i < test_size; i++) {
        auto item = umap[uniform_dist(e1)];
        sum += item.id_;
    }
}
*/
BENCHMARK_MAIN();
